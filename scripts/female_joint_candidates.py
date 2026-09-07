#!/usr/bin/env python3
"""Fit disabled HRA knee/tendon candidates against retained atlas joint surfaces.

Requires numpy. Writes only --output-dir. Surface crops/proximity patches are
algorithmic regions, NOT reviewed anatomical landmarks or validated entheses.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('joint_surface_helpers', ROOT/'scripts/pelvis-surface-audit.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
BONES = {'L': {'femur':'FJ3259','tibia':'FJ3282','patella':'FJ3275'},
         'R': {'femur':'FJ3365','tibia':'FJ3387','patella':'FJ3381'}}
STRUCTURES = ['anterior_cruciate_ligament_of_knee', 'posterior_cruciate_ligament_of_knee',
              'anterolateral_ligament_of_knee', 'meniscus', 'tendon_of_quadriceps_femoris']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs(root):
    paths = {root/'public/models/female-fit-report.json'}
    for name in ['atlas.json','atlas-female.json','atlas-female-reconstructed.json']:
        path = root/'public/models'/name
        paths.add(path)
        for c in json.loads(path.read_text())['chunks']:
            paths.add(root/'public'/c['url'].lstrip('/'))
    return {str(p.relative_to(root)):sha(p) for p in sorted(paths)}


def sample(points, count=72, phase=0):
    # Stable disjoint even/odd pools for fitting and holdout where large enough.
    pool = np.arange(phase, len(points), 2)
    return points[pool[np.linspace(0,len(pool)-1,min(count,len(pool)),dtype=int)]]


def crop(v, f, name):
    # Only triangles fully inside the geometric local-height region are retained.
    mask = np.ones(len(v),dtype=bool)
    if name == 'femur': mask = v[:,1] <= v[:,1].min()+.09
    if name == 'tibia': mask = v[:,1] >= v[:,1].max()-.08
    faces = f[mask[f].all(axis=1)]
    used = np.unique(faces)
    assert len(faces)>0
    return v[used], a.Surface(v,faces)


def similarity(source, target):
    x0,y0 = source.mean(axis=0),target.mean(axis=0)
    x,y=source-x0,target-y0
    u,_,vt=np.linalg.svd(x.T@y)
    correction=np.eye(3); correction[-1,-1]=np.linalg.det(u@vt)
    rotation=u@correction@vt  # row vectors
    scale=float(np.clip(np.sum((x@rotation)*y)/np.sum(x*x), .85, 1.20))
    return scale, rotation, y0-scale*x0@rotation


def fit_joint(source_meshes,target_meshes):
    source_regions={k:crop(*v,k) for k,v in source_meshes.items()}
    target_regions={k:crop(*v,k) for k,v in target_meshes.items()}
    cent_s=np.array([source_regions[k][0].mean(axis=0) for k in BONES['L']])
    cent_t=np.array([target_regions[k][0].mean(axis=0) for k in BONES['L']])
    # Translation initialization retains already verified source axis orientation.
    scale,rotation,translation=1.,np.eye(3),(cent_t-cent_s).mean(axis=0)
    fitting={k:sample(v[0]) for k,v in source_regions.items()}
    history=[]
    for iteration in range(18):
        src=[]; dst=[]; distances=[]
        for k,points in fitting.items():
            q=scale*points@rotation+translation
            nearest=[target_regions[k][1].nearest(p) for p in q]
            src.append(points);dst.append(np.array([x[1] for x in nearest]))
            distances.extend(x[0] for x in nearest)
        history.append(float(np.sqrt(np.mean(np.square(distances)))*1000))
        scale,rotation,translation=similarity(np.concatenate(src),np.concatenate(dst))
    transform=lambda points:scale*points@rotation+translation
    heldout={}
    for k,(points,_) in source_regions.items():
        heldout[k]=a.stats(target_regions[k][1].distances(transform(sample(points,96,1))))
    return transform,dict(scale=scale,rotationRowVectors=rotation.tolist(),translationM=translation.tolist(),
                          determinant=float(np.linalg.det(scale*rotation)),iterations=18,
                          fittingRmsMm=history,heldOutSourceToTargetSurface=heldout,
                          sourceCropVertexCounts={k:len(v[0]) for k,v in source_regions.items()})


def refine_joint(source_meshes,target_meshes,baseline,baseline_fit,ridge=.015):
    """One Gaussian RBF displacement field constrained by local bone ICP fits.

    Correspondences are automatic geometric surface points, not annotated anatomy.
    Candidate soft-tissue vertices are NOT used to fit this field.
    """
    controls=[]; corrections=[]; local={}
    for bone in source_meshes:
        source_points,_=crop(*source_meshes[bone],bone)
        _,target_surface=crop(*target_meshes[bone],bone)
        points=sample(source_points,72)
        scale=baseline_fit['scale']; rotation=np.array(baseline_fit['rotationRowVectors'])
        translation=np.array(baseline_fit['translationM'])
        for _ in range(16):
            current=scale*points@rotation+translation
            nearest=np.array([target_surface.nearest(p)[1] for p in current])
            scale,rotation,translation=similarity(points,nearest)
        current=scale*points@rotation+translation
        nearest=np.array([target_surface.nearest(p)[1] for p in current])
        # The local similarity preserves correspondence orientation; nearest surface
        # correction addresses donor-shape mismatch. Regularization smooths it.
        controls.append(points);corrections.append(nearest-baseline(points))
        local[bone]=dict(scale=scale,rotationRowVectors=rotation.tolist(),translationM=translation.tolist(),
                         localSimilaritySampleResidual=a.stats(target_surface.distances(current)))
    centers=np.concatenate(controls); residual=np.concatenate(corrections)
    radius=.025
    kernel=np.exp(-np.sum((centers[:,None]-centers[None,:])**2,axis=2)/(2*radius**2))
    coefficients=np.linalg.solve(kernel+ridge*np.eye(len(centers)),residual)
    def field(points):
        kernel=np.exp(-np.sum((points[:,None]-centers[None,:])**2,axis=2)/(2*radius**2))
        return baseline(points)+kernel@coefficients
    def jacobian(points):
        offset=points[:,None]-centers[None,:]
        weights=np.exp(-np.sum(offset*offset,axis=2)/(2*radius**2))
        # J rows output components, columns input coordinates.
        return baseline_fit['scale']*np.array(baseline_fit['rotationRowVectors']).T[None,:,:]-np.einsum('nc,ncj,ci->nij',weights,offset,coefficients)/(radius**2)
    holdout={}
    for bone,mesh in source_meshes.items():
        points,_=crop(*mesh,bone); _,surface=crop(*target_meshes[bone],bone)
        holdout[bone]=a.stats(surface.distances(field(sample(points,96,1))))
    metadata=dict(method='Bone-local similarity ICP then one regularized Gaussian RBF surface displacement',
        localBoneFits=local,kernelRadiusM=radius,ridge=ridge,controlCount=len(centers),
        controlSourcePointsM=centers.tolist(),controlDisplacementsM=residual.tolist(),coefficients=coefficients.tolist(),
        heldOutSourceToTargetSurface=holdout,
        limitations='Automatic nearest-surface correspondences; no reviewed landmark or endpoint annotations; no soft-tissue points used in fitting.')
    return field,jacobian,metadata


def distortion(v,f,q,baseline_q,jacobian):
    edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
    length=np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1)
    qlength=np.linalg.norm(q[edges[:,1]]-q[edges[:,0]],axis=1)
    blength=np.linalg.norm(baseline_q[edges[:,1]]-baseline_q[edges[:,0]],axis=1)
    mask=(length>1e-9)&(blength>1e-9)
    points=np.concatenate([v,v[f].mean(axis=1)])
    jac=jacobian(points)
    det=np.linalg.det(jac); singular=np.linalg.svd(jac,compute_uv=False)
    return dict(evaluatedVerticesAndTriangleCentroids=len(points),minimumJacobianDeterminant=float(det.min()),
        nonpositiveJacobianSamples=int((det<=0).sum()),minimumSingularValue=float(singular.min()),maximumSingularValue=float(singular.max()),
        finalEdgeLengthRatioToSource=np.quantile(qlength[mask]/length[mask],[0,.05,.5,.95,1]).tolist(),
        finalEdgeLengthRatioToSimilarity=np.quantile(qlength[mask]/blength[mask],[0,.05,.5,.95,1]).tolist(),
        note='Jacobian is for pre-morph refinement, edges include unchanged female morph; sampled positivity is not a global injectivity proof.')


def connected_components(v,f):
    _,indices=np.unique(np.round(v,7),axis=0,return_inverse=True)
    parents=np.arange(indices.max()+1)
    def find(x):
        while parents[x]!=x:
            parents[x]=parents[parents[x]];x=parents[x]
        return x
    for tri in indices[f]:
        for j in [1,2]:parents[find(tri[j])]=find(tri[0])
    return len({find(int(x)) for x in indices})


def export_mesh(output,key,v,f):
    normal=np.zeros_like(v)
    tri=v[f]; face_normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    for j in range(3):np.add.at(normal,f[:,j],face_normal)
    norm=np.linalg.norm(normal,axis=1);normal/=np.maximum(norm[:,None],1e-30)
    path=output/f'{key}.npz'
    # .npz container timestamps are not a reproducibility contract; raw arrays are.
    np.savez(path,positions=v.astype('<f4'),indices=f.astype('<u4'),normals=normal.astype('<f4'))
    obj=output/f'{key}.obj'
    with obj.open('w') as handle:
        for p in v:handle.write('v %.9f %.9f %.9f\n'%tuple(p))
        for t in f:handle.write('f %d %d %d\n'%tuple(t+1))
    return dict(npz=path.name,obj=obj.name,vertices=len(v),triangles=len(f),
                positionsSha256=hashlib.sha256(v.astype('<f4').tobytes()).hexdigest(),
                indicesSha256=hashlib.sha256(f.astype('<u4').tobytes()).hexdigest(),
                boundsM=[v.min(axis=0).tolist(),v.max(axis=0).tolist()])


def diagnostic_plot(root,output):
    """Orthographic triangle projections, not a rendered anatomical approval."""
    female=a.Atlas(root,'atlas-female-reconstructed.json')
    colors=['#df554c','#486bc6','#db8b18','#289b86','#a855c7']
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="1250" viewBox="0 0 1400 1250">',
         '<rect width="1400" height="1250" fill="white"/>',
         '<text x="30" y="30" font-size="22">HRA joint candidates — UNAPPROVED geometric fit</text>',
         '<text x="30" y="55" font-size="14">Retained bones gray; exact mesh triangle projections. Cropped knee region; anterior is right in sagittal views.</text>']
    for side_index,side in enumerate(['L','R']):
        bones={k:female.mesh([v]) for k,v in BONES[side].items()}
        center=bones['patella'][0].mean(axis=0)
        for view_index,(name,axis) in enumerate([('Frontal',0),('Sagittal',2)]):
            ox=side_index*700;oy=90+view_index*545
            svg.append(f'<text x="{ox+20}" y="{oy+20}" font-size="19">{side} — {name}</text>')
            # Same physical scale in all panels, 2400 pixels per meter.
            def polygons(vertices,faces,color,opacity):
                triangles=vertices[faces]
                keep=(triangles[:,:,1].max(axis=1)>center[1]-.10)&(triangles[:,:,1].min(axis=1)<center[1]+.105)
                triangles=triangles[keep]
                for tri in triangles:
                    xy=np.column_stack([ox+350+(tri[:,axis]-center[axis])*2400,oy+280-(tri[:,1]-center[1])*2400])
                    points=' '.join(f'{x:.2f},{y:.2f}' for x,y in xy)
                    svg.append(f'<polygon points="{points}" fill="{color}" fill-opacity="{opacity}" stroke="{color}" stroke-opacity=".25" stroke-width=".35"/>')
            svg.append(f'<defs><clipPath id="clip{side}{view_index}"><rect x="{ox+10}" y="{oy+28}" width="680" height="510"/></clipPath></defs>')
            svg.append(f'<g clip-path="url(#clip{side}{view_index})">')
            for v,f in bones.values():polygons(v,f,'#7c8891',.045)
            for structure,color in zip(STRUCTURES,colors):
                mesh=np.load(output/f'VH_F_{structure}_{side}.npz')
                polygons(mesh['positions'],mesh['indices'],color,.20)
            svg.append('</g>')
    for i,(name,color) in enumerate(zip(['ACL','PCL','Anterolateral ligament','Menisci','Quadriceps tendon'],colors)):
        x=30+i*270
        svg.append(f'<rect x="{x}" y="1210" width="18" height="18" fill="{color}"/><text x="{x+25}" y="1225" font-size="14">{name}</text>')
    svg.append('</svg>')
    (output/'joint-diagnostic.svg').write_text('\n'.join(svg))


def run(root,output,regularization=.015):
    if not np.isfinite(regularization) or regularization<=0:raise ValueError("Regularization must be finite and positive")
    root=root.resolve();output=output.resolve()
    if output==ROOT.resolve() or ROOT.resolve() in output.parents or output==root or root in output.parents:
        raise ValueError('Candidate output must resolve outside the repository, including symlinks')
    output.mkdir(parents=True,exist_ok=True)
    before=inputs(root)
    source=a.Atlas(root,'atlas-female.json'); target=a.Atlas(root,'atlas.json')
    female=a.Atlas(root,'atlas-female-reconstructed.json')
    morph_settings=json.loads((root/'public/models/female-fit-report.json').read_text())['morph']
    result=dict(schemaVersion=1,enabled=False,status='candidate-only; attachment annotations unreviewed',
                inputHashes=before,method={'registration':'Shared similarity baseline followed by per-bone ICP and one regularized Gaussian RBF displacement field',
                'sourceFrame':'Bundled Y-up meters; +Z anterior; no axis flip',
                'surfaceCrops':'Femur includes its 15 source compound leaves excluding cartilage. All patella triangles; femur bottom 90 mm; tibia top 80 mm; crops are not anatomical landmarks',
                'scaleBounds':[.85,1.2],'attachmentPatchThresholdMm':3,
                'sampling':'Sampled surface fitting and holdout only, not complete collision validation. 72 even-index source region vertices/bone fit; up to96 odd-index holdout; deterministic evenly spaced indices',
                'morph':'Apply existing fit-report morph once after registration; body is unchanged',
                'normalMethod':'Area-weighted normals recomputed from candidate triangles',
                'limitations':['One-way surface ICP can have local minima; fitting regions and attachment patches are unreviewed.',
                'Unsigned distances do not establish tissue containment, intersection absence, or functional attachment.',
                'HRA knees are approximately mirrored assembly parts; this is not a validated female-specific joint fit.']},
                candidates={'knee':{'enabled':False,'status':'blocked: joint fit and attachment placement require refinement','partIds':[]},'quadricepsTendons':{'enabled':False,'status':'blocked: patellar interface does not fit','partIds':[]}},sides={})
    concepts={c['id']:c for c in source.data['concepts']}
    for side in ['L','R']:
        print(f'Fitting {side} joint',flush=True)
        source_ids={k:[f'VH_F_{k}_{side}'] for k in BONES[side]}
        # The main femur leaf omits separately segmented condylar/enthesis regions.
        # Use the source-defined compound, explicitly excluding its cartilage leaf.
        source_ids['femur']=sorted(key for key in concepts[f'HRA:VH_F_femur_{side}']['elements'] if 'cartilage' not in key)
        assert len(source_ids['femur'])==15
        source_bones={k:source.mesh(ids) for k,ids in source_ids.items()}
        target_bones={k:target.mesh([key]) for k,key in BONES[side].items()}
        baseline,fit=fit_joint(source_bones,target_bones)
        transform,jacobian,refinement=refine_joint(source_bones,target_bones,baseline,fit,regularization)
        source_surfaces={k:a.Surface(*v) for k,v in source_bones.items()}
        female_surfaces={k:a.Surface(*female.mesh([key])) for k,key in BONES[side].items()}
        rows=[]
        for structure in STRUCTURES:
            key=f'VH_F_{structure}_{side}'
            group='quadricepsTendons' if structure.startswith('tendon') else 'knee'
            result['candidates'][group]['partIds'].append(key)
            v,f=source.mesh([key]);q=a.morph(transform(v),morph_settings,key)
            baseline_q=a.morph(baseline(v),morph_settings,key)
            row=dict(id=key,name=source.parts[key]['name'],candidate=group,export=export_mesh(output,key,q,f),
                     topologyPreserved=True,connectedComponents=connected_components(v,f),boneScreens={},
                     distortion=distortion(v,f,q,baseline_q,jacobian))
            for bone,surface in source_surfaces.items():
                sd=surface.distances(v); mask=sd<=.003
                fd=female_surfaces[bone].distances(q)
                bd=female_surfaces[bone].distances(baseline_q[mask]) if mask.any() else None
                row['boneScreens'][bone]=dict(sourceAllVertices=a.stats(sd),candidateAllVertices=a.stats(fd),
                    sourceProximityPatchVertexCount=int(mask.sum()),
                    sourcePatch=a.stats(sd[mask]) if mask.any() else None,
                    candidateSamePatch=a.stats(fd[mask]) if mask.any() else None,
                    similaritySamePatch=a.stats(bd) if bd is not None else None,
                    note='Source vertices within 3 mm of this bone; not reviewed insertion/enthesis vertices')
            if group=='quadricepsTendons':
                src_muscle=a.Surface(*source.mesh([f'VH_F_rectus_femoris_{side}']))
                sd=src_muscle.distances(v);mask=sd<=.003
                muscle_ids=[p+('M' if side=='L' else '') for p in ['FJ1433','FJ1441','FJ1442','FJ1443']]
                td=a.Surface(*female.mesh(muscle_ids)).distances(q)
                row['quadricepsInterface']=dict(targetIds=muscle_ids,sourceRectusAllVertices=a.stats(sd),
                    candidateQuadricepsAllVertices=a.stats(td),sourceNearRectusCount=int(mask.sum()),
                    candidateSamePatch=a.stats(td[mask]) if mask.any() else None,
                    similaritySamePatch=a.stats(a.Surface(*female.mesh(muscle_ids)).distances(baseline_q[mask])) if mask.any() else None,
                    candidateVerticesWithin1Mm=int((td<=.001).sum()),candidateVerticesWithin3Mm=int((td<=.003).sum()),
                    duplicationConclusion='Proximity alone cannot distinguish intended joining from duplicate embedded distal tendon; visual/semantic review required')
            rows.append(row)
            print(f'  screened {key}',flush=True)
        result['sides'][side]=dict(sourceBoneRegionIds=source_ids,targetBoneIds=BONES[side],fit=fit,refinement=refinement,parts=rows)
    after=inputs(root)
    assert before==after,'Input atlas changed during run; reject this report and rerun.'
    result['existingAtlasByteIdentityVerified']=True
    result['candidates']['knee']['blockers']=['Target attachment landmarks are not reviewed; residuals are geometric proximity evidence only.',
        'Meniscus two-component labels and ligament footprint placement require joint section review before enabling.']
    result['candidates']['quadricepsTendons']['blockers']=['Depends on accepted knee fit.',
        'Existing target quadriceps may already contain distal tendinous geometry; duplication and four-head junction remain unreviewed.']
    diagnostic_plot(root,output)
    result['diagnosticPlot']='joint-diagnostic.svg'
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def self_test():
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    angle=.2;r=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    expected=points@r*1.1+np.array([.1,.2,.3])
    s,rr,t=similarity(points,expected)
    assert np.allclose(s*points@rr+t,expected,atol=1e-12)
    assert np.linalg.det(rr)>0
    surface=a.Surface(np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]),np.array([[0,1,2]]))
    assert np.allclose(surface.nearest(np.array([.2,.2,1.]))[0],1.)
    assert np.allclose(surface.nearest(np.array([2.,0.,0.]))[0],1.)
    cube=np.array([[x,y,z] for x in [0.,.01] for y in [0.,.01] for z in [0.,.01]])
    faces=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
    meshes={bone:(cube+np.array([i*.025,0,0]),faces) for i,bone in enumerate(['femur','tibia','patella'])}
    targets={bone:(v+np.array([0,.001,.002]),f) for bone,(v,f) in meshes.items()}
    field,jac,_=refine_joint(meshes,targets,lambda x:x,dict(scale=1.,rotationRowVectors=np.eye(3).tolist(),translationM=[0,0,0]))
    probes=np.array([[.006,.004,.005],[.03,.013,.008]])
    numerical=np.stack([(field(probes+np.eye(3)[axis]*1e-6)-field(probes-np.eye(3)[axis]*1e-6))/(2e-6) for axis in range(3)],axis=2)
    assert np.allclose(jac(probes),numerical,atol=1e-7)
    assert np.array_equal(field(np.repeat(probes[:1],2,axis=0))[0],field(probes[:1])[0])
    assert np.all(np.linalg.det(jac(probes))>0)
    print('Similarity recovery, exact face/edge distances, coherent RBF field and analytical Jacobian checks pass.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output-dir',type=Path,default=Path('/tmp/female-joint-candidates'))
    parser.add_argument('--regularization',type=float,default=.015,help='Gaussian RBF ridge; .015 historical refined candidate, .08 smoother ALL preview')
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--plot-only',action='store_true',help='Plot an existing report after checking current input hashes')
    args=parser.parse_args()
    if args.self_test:self_test()
    elif args.plot_only:
        resolved=args.output_dir.resolve()
        for protected in [ROOT.resolve(),args.root.resolve()]:
            if resolved==protected or protected in resolved.parents:
                raise ValueError('Candidate output must resolve outside the repository, including symlinks')
        report=json.loads((args.output_dir/'report.json').read_text())
        assert report['inputHashes']==inputs(args.root)
        diagnostic_plot(args.root,args.output_dir)
    else:run(args.root,args.output_dir,args.regularization)
