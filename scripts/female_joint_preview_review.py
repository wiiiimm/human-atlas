#!/usr/bin/env python3
"""Independent geometry checks for the explicitly selected bilateral ALL preview.

Writes a separate report into the isolated candidate directory. Does not modify
candidate meshes or the candidate fit report. Requires numpy only.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.dont_write_bytecode=True
import female_joint_candidates as j

IDS=[f'VH_F_anterolateral_ligament_of_knee_{s}' for s in ['L','R']]


def closed(v,f):
    _,mapping=np.unique(np.round(v,7),axis=0,return_inverse=True)
    faces=mapping[f]
    edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    _,counts=np.unique(edges,axis=0,return_counts=True)
    return dict(boundaryEdges=int((counts==1).sum()),nonmanifoldEdges=int((counts>2).sum()),weldPrecisionM=1e-7)


def classify(point,triangles):
    aa,bb,cc=triangles[:,0],triangles[:,1],triangles[:,2]
    e1,e2=bb-aa,cc-aa;votes=[]
    for direction in [np.array([1.,.331,.127]),np.array([.173,1.,.417]),np.array([.231,.377,1.])]:
        h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h);valid=abs(det)>1e-14
        inv=np.divide(1,det,out=np.zeros_like(det),where=valid);s=point-aa
        u=inv*np.einsum('ij,ij->i',s,h);q=np.cross(s,e1)
        v=inv*(q@direction);t=inv*np.einsum('ij,ij->i',e2,q)
        hit=valid&(u>=0)&(v>=0)&(u+v<=1)&(t>1e-9)
        if np.any(hit&((u<1e-8)|(v<1e-8)|(1-u-v<1e-8))):continue
        hits=np.sort(t[hit]);unique=hits[np.r_[True,np.diff(hits)>1e-7]] if len(hits) else hits
        votes.append(len(unique)%2)
    return votes[0] if len(votes)>=2 and len(set(votes))==1 else -1


def segment_hits(p0,p1,tri):
    e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];direction=p1-p0
    h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h);valid=abs(det)>1e-16
    inv=np.divide(1,det,out=np.zeros_like(det),where=valid);s=p0-tri[:,0]
    u=inv*np.einsum('ij,ij->i',s,h);q=np.cross(s,e1)
    v=inv*np.einsum('ij,ij->i',direction,q);t=inv*np.einsum('ij,ij->i',e2,q)
    return valid&(u>=0)&(v>=0)&(u+v<=1)&(t>1e-8)&(t<1-1e-8)


def crossings(v,f,bv,bf):
    target=bv[bf];lo=target.min(axis=1);hi=target.max(axis=1);count=0;pairs=0
    for tri in v[f]:
        select=((hi>=tri.min(axis=0))&(lo<=tri.max(axis=0))).all(axis=1)
        other=target[select]
        if not len(other):continue
        hits=np.zeros(len(other),dtype=bool);repeated=np.broadcast_to(tri,other.shape)
        for i in range(3):
            hits|=segment_hits(np.broadcast_to(tri[i],(len(other),3)),np.broadcast_to(tri[(i+1)%3],(len(other),3)),other)
            hits|=segment_hits(other[:,i],other[:,(i+1)%3],repeated)
        count+=int(hits.any());pairs+=int(hits.sum())
    return dict(tissueTrianglesWithNoncoplanarSurfaceCrossing=count,intersectingTrianglePairs=pairs,
                method='AABB broad phase; both triangles edge/face tests; coplanar contacts and exact endpoint-only contacts excluded')


def intrusion(v,f,bv,bf):
    topology=closed(bv,bf)
    result={'boneClosure':topology,'triangleCrossings':crossings(v,f,bv,bf)}
    if topology['boundaryEdges'] or topology['nonmanifoldEdges']:
        result['insideClassification']='unavailable: bone is not a closed two-manifold after position welding'
        return result
    points=np.concatenate([v,v[f].mean(axis=1)])
    flags=np.array([classify(p,bv[bf]) for p in points]);inside=points[flags==1]
    depths=j.a.Surface(bv,bf).distances(inside) if len(inside) else np.array([])
    result.update(sampleCount=len(points),insideCount=len(inside),ambiguousCount=int((flags==-1).sum()),
        insideDepth=j.a.stats(depths) if len(depths) else None,
        note='Every vertex and triangle centroid screened with three ray directions; depth is exact point-to-triangle distance, not continuous maximum penetration.')
    return result


def run(root,candidate):
    root=root.resolve();candidate=candidate.resolve()
    if candidate==root or root in candidate.parents or candidate==j.ROOT or j.ROOT in candidate.parents:
        raise ValueError('Candidate must resolve outside repository')
    report_path=candidate/'report.json';fit=json.loads(report_path.read_text())
    assert fit['inputHashes']==j.inputs(root)
    source=j.a.Atlas(root,'atlas-female.json');female=j.a.Atlas(root,'atlas-female-reconstructed.json')
    result=dict(schemaVersion=1,scope='Experimental static study preview: bilateral anterolateral ligament only',
        recommendedIds=IDS,candidateReportSha256=j.sha(report_path),inputHashes=fit['inputHashes'],parts=[])
    for side in ['L','R']:
        key=f'VH_F_anterolateral_ligament_of_knee_{side}';data=fit['sides'][side]
        row=next(p for p in data['parts'] if p['id']==key)
        path=candidate/row['export']['npz'];mesh=np.load(path);v=mesh['positions'].astype(float);f=mesh['indices'];n=mesh['normals']
        sv,sf=source.mesh([key]);assert np.array_equal(sf,f)
        assert np.isfinite(v).all() and np.isfinite(n).all()
        assert np.allclose(np.linalg.norm(n,axis=1),1,atol=1e-5)
        assert np.allclose([v.min(axis=0),v.max(axis=0)],row['export']['boundsM'],atol=1e-7,rtol=0)
        assert j.hashlib.sha256(mesh['positions'].astype('<f4').tobytes()).hexdigest()==row['export']['positionsSha256']
        assert j.hashlib.sha256(f.astype('<u4').tobytes()).hexdigest()==row['export']['indicesSha256']
        expected=np.zeros_like(v);tri=v[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        for axis in range(3):np.add.at(expected,f[:,axis],fn)
        expected/=np.maximum(np.linalg.norm(expected,axis=1)[:,None],1e-30)
        # Export normals computed before float32 quantization; tolerate tiny differences.
        assert np.allclose(expected,n,atol=5e-4)
        part=dict(id=key,sourceTopologyExact=True,finitePositionsAndNormals=True,unitNormals=True,
            normalsMatchDeformedSurface=True,boundsMatchReport=True,npzSha256=j.sha(path),
            normalizedNormalsSha256=j.hashlib.sha256(n.astype('<f4').tobytes()).hexdigest(),
            shapeDistortion=row['distortion'],sourceVsCandidateAttachmentDistances=row['boneScreens'],boneIntrusions={})
        for bone in ['femur','tibia']:
            bv,bf=female.mesh([data['targetBoneIds'][bone]])
            sbv,sbf=source.mesh(data['sourceBoneRegionIds'][bone])
            part['boneIntrusions'][bone]=dict(candidate=intrusion(v,f,bv,bf),source=intrusion(sv,sf,sbv,sbf))
        result['parts'].append(part)
    assert fit['inputHashes']==j.inputs(root)
    result['existingAtlasUnchanged']=True
    (candidate/'all-preview-review.json').write_text(json.dumps(result,indent=2)+'\n')
    print(candidate/'all-preview-review.json')


def self_test():
    v=np.array([[x,y,z] for x in [0.,1.] for y in [0.,1.] for z in [0.,1.]])
    f=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
    assert closed(v,f)['boundaryEdges']==0
    assert classify(np.array([.5,.5,.5]),v[f])==1
    assert classify(np.array([2.,.5,.5]),v[f])==0
    tri=np.array([[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]])
    assert segment_hits(np.array([[.2,.2,-1.]]),np.array([[.2,.2,1.]]),tri)[0]
    assert not segment_hits(np.array([[2.,2.,-1.]]),np.array([[2.,2.,1.]]),tri)[0]
    print('Closed-cube parity and triangle segment hit/miss tests passed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=j.ROOT)
    parser.add_argument('--candidate-dir',type=Path,default=Path('/tmp/female-joint-preview-candidates'))
    parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:self_test()
    else:run(args.root,args.candidate_dir)
