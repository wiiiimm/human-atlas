"""Build a female study prototype from the BodyParts3D framework and HRA female organs.

Run: python3 scripts/build-female-reconstruction.py (requires numpy).

Pipeline
1. Omit male-specific BodyParts3D structures from the manifest.
2. Fit HRA reproductive organs and female pelvis into the retained framework.
3. Apply the recorded body morph to non-breast structures, including explicitly
   scoped glute contours. A finite Jacobian sample is not anatomical validation.
4. Transform all 16 original HRA breast meshes directly into the final atlas frame:
   retain the base footprint, reduce anterior projection, and lift lower-front tissue.
   Breast meshes do not receive the body morph or the former drape/inset machinery.

5. Append the recorded bilateral anterolateral knee ligament fit, retaining all
   existing model buffers. The candidate recipe is pinned in data/anatomy/joint-additions.

Parameters and source hashes are recorded in the fit report. Source topology is
preserved; placement and proportions remain experimental and unreviewed anatomy.
"""
from pathlib import Path
import argparse, copy, glob, gzip, hashlib, json, os, re, subprocess, shutil, sys, tempfile
import numpy as np
import hra_breast_contour as breast_contour

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'public/models'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=INPUT,help='Write a candidate to a separate directory before publishing it.')
OUT=parser.parse_args().output_dir.resolve()
OUT.mkdir(parents=True,exist_ok=True)
male=json.loads((INPUT/'atlas.json').read_text())
female=json.loads((INPUT/'atlas-female.json').read_text())
mp={p['id']:p for p in male['parts']}
fp={p['id']:p for p in female['parts']}
male_buffers=[(ROOT/('public'+c['url'])).read_bytes() for c in male['chunks']]
source_buffers=[(ROOT/('public'+c['url'])).read_bytes() for c in female['chunks']]

def read(buffers,p):
    b=buffers[p['chunk']]
    pos=np.frombuffer(b,'<f4',p['vertexCount']*3,p['positions']).reshape(-1,3).astype(np.float64)
    nrm=np.frombuffer(b,'<i2',p['vertexCount']*3,p['normals']).reshape(-1,3).astype(np.float64)/32767
    idx=np.frombuffer(b,'<u4',p['indexCount'],p['indices']).copy()
    return pos,nrm,idx

# ---------------------------------------------------------------- 1. exclusions
excluded_pattern=re.compile(r'penis|testicular|testis|scrot|prostat|seminal|deferent|epididym|spermatic|ejaculat|cremaster|bulbospong|ischiocav|perineal|coccygeus|puborectalis|levator ani|sphincter',re.I)
# The male pelvis is replaced by the HRA female pelvis (compact bone shells, sacrum, coccyx).
replacements={
 'FJ3152':['VH_F_ilium_compact_bone_R','VH_F_ischium_compact_bone_R','VH_F_pubis_compact_bone_R'],
 'FJ3288':['VH_F_ilium_compact_bone_L','VH_F_ischium_compact_bone_L','VH_F_pubis_compact_bone_L'],
 'FJ3393':['VH_F_sacrum','VH_F_coccyx'],
}
pelvis_names={'VH_F_ilium_compact_bone_R':'Right ilium','VH_F_ischium_compact_bone_R':'Right ischium','VH_F_pubis_compact_bone_R':'Right pubis','VH_F_ilium_compact_bone_L':'Left ilium','VH_F_ischium_compact_bone_L':'Left ischium','VH_F_pubis_compact_bone_L':'Left pubis','VH_F_sacrum':'Sacrum','VH_F_coccyx':'Coccyx'}
retained=[];excluded=[]
for p in male['parts']:
    reason=('Replaced by the HRA female pelvis' if p['id'] in replacements else
            'Male body surface omitted' if p['system']=='integumentary' else
            'Male reproductive anatomy omitted' if p['system']=='reproductive' else
            'Male urethra omitted; no validated female replacement bundled' if p['id']=='FJ3148' else
            'Sex-specific or pelvic-floor structure omitted' if excluded_pattern.search(p['name']) else None)
    if reason:excluded.append(dict(id=p['id'],name=p['name'],reason=reason))
    else:retained.append(copy.deepcopy(p))

def smoothstep(a,b,t):
    u=np.clip((t-a)/(b-a),0,1);return u*u*(3-2*u)
def bounds_of(parts):
    b=np.array([p['bounds'] for p in parts]);return b[:,0].min(axis=0),b[:,1].max(axis=0)
def bounds_fit(source_parts,target_parts):
    lo,hi=bounds_of(source_parts);tlo,thi=bounds_of(target_parts)
    scale=(thi-tlo)/(hi-lo);return scale,tlo-lo*scale

# ---------------------------------------------------------------- 2. affine fits
# One transform for the whole reproductive assembly preserves its internal alignment. The HRA
# bladder envelope is matched to the retained bladder as a placement proxy.
bladder_ids=['VH_F_fundus_of_urinary_bladder_dome','VH_F_fundus_of_urinary_bladder_base']
bladder_bounds=np.array([fp[i]['bounds'] for i in bladder_ids])
lo=bladder_bounds[:,0].min(axis=0);hi=bladder_bounds[:,1].max(axis=0)
target=np.array(mp['FJ3149']['bounds'])
pelvic_scale=(target[1]-target[0])/(hi-lo)
pelvic_offset=target[0]-lo*pelvic_scale
# Female pelvis: match hip bone envelopes as an initial placement proxy.
# This does not establish acetabular contact, landmark alignment, or muscle attachments.
hip_ids=[i for k in ('FJ3152','FJ3288') for i in replacements[k]]
pelvis_scale,pelvis_offset=bounds_fit([fp[i] for i in hip_ids],[mp['FJ3152'],mp['FJ3288']])
transforms={
 'reproductive':dict(scale=pelvic_scale.tolist(),offset=pelvic_offset.tolist(),basis='HRA bladder bounds aligned to retained BodyParts3D bladder bounds'),
 'pelvis':dict(scale=pelvis_scale.tolist(),offset=pelvis_offset.tolist(),basis='HRA hip bone bounds aligned to BodyParts3D hip bone bounds; sacrum and coccyx share the transform')
}
pelvis_ids=[i for ids in replacements.values() for i in ids]
selected=[p for p in female['parts'] if p['system']=='reproductive' or (p['system']=='integumentary' and p['id']!='VH_F_skin') or p['id'] in pelvis_ids]

# Assemble every mesh in the base coordinate frame before deformation.
meshes=[]  # dict(part, pos, nrm, idx, provenance)
for p in retained:
    pos,nrm,idx=read(male_buffers,p)
    meshes.append(dict(part=p,pos=pos,nrm=nrm,idx=idx,provenance={'source':'BodyParts3D 4.0','sourceId':p['id'],'adaptation':'Retained reference geometry; reshaped by the shared female body morph'}))
contoured=[]
for source in selected:
    group='pelvis' if source['id'] in pelvis_ids else 'mammary' if source['id'] in breast_contour.BREAST_IDS else 'reproductive'
    pos,nrm,idx=read(source_buffers,source)
    if group!='mammary':
        t=transforms[group];scale=np.array(t['scale']);offset=np.array(t['offset'])
        pos=pos*scale+offset
        nrm=nrm/scale;nrm/=np.maximum(np.linalg.norm(nrm,axis=1,keepdims=True),1e-20)
    p=copy.deepcopy(source);p['system']='skeletal' if group=='pelvis' else group
    if group=='pelvis':p['name']=pelvis_names[source['id']]
    if group=='mammary':
        contoured.append(source['id'])
        if source['id'] in ('VH_F_fat_L','VH_F_fat_R'):
            p['name']='Adipose tissue of '+('left' if source['id'].endswith('_L') else 'right')+' breast'
        adaptation='Original HRA breast mesh; source base footprint retained, anterior projection reduced and lower front lifted; translated directly into final female coordinates without body morph. Experimental placement.'
    else:
        adaptation=('Female pelvis fitted to the BodyParts3D hip bone envelope, replacing the male pelvis' if group=='pelvis' else 'Female reference anatomy fitted to the BodyParts3D framework')+'; experimental placement reshaped by the shared female body morph'
    meshes.append(dict(part=p,pos=pos,nrm=nrm,idx=idx,group=group,provenance={'source':'HRA united-female v1.5','sourceId':source['id'],'adaptation':adaptation}))
assert set(contoured)==set(breast_contour.BREAST_IDS)

# ---------------------------------------------------------------- 3. non-breast body morph
# Lateral width factor by height (base frame, before stature scaling). Between knots the
# factor follows a smoothstep so the field is continuous and slopes stay gentle along the arms.
MORPH=dict(
 stature=.95,
 # User-reference silhouette pass: broader pelvis/upper hips and narrower shoulder cap.
 # Artistic proportions, not a population measurement. The Sketchfab/comparison waist pass
 # keeps the slim waist; a three-view pass rebalances hip width and posterior glute contour.
 lateralKnots=[[0.,1.],[.55,1.],[.80,1.08],[.92,1.15],[1.03,1.15],[1.12,.955],[1.21,.96],[1.29,.93],[1.42,.88],[1.48,.88],[1.60,1.]],
 # Retain the existing limb field outside the torso, with a nearly rigid 3 mm inward arm shift.
 # Blend the new torso silhouette smoothly rather than sending hip widening into the forearms.
 lateralBlend=dict(referenceKnots=[[0.,1.],[.55,1.],[.80,1.05],[.92,1.10],[1.03,1.10],[1.21,.96],[1.29,.95],[1.42,.91],[1.48,.91],[1.60,1.]],innerRadius=.17,outerRadius=.23,outerTranslation=-.003,heightRamp=[.45,.55,1.48,1.60]),
 waistRefinement=dict(delta=-.095,heightRamp=[1.03,1.12,1.21],radialRamp=[.14,.17]),
 lateralRadius=.17,  # beyond this |x| the lateral change saturates into a translation, so arms move with the torso instead of being squeezed
 thoraxDepth=dict(scale=.96,ramp=[1.05,1.15,1.35,1.45]),
 # Posterior and inferior contour for explicitly listed glute muscles and their inferior veins.
 # Paired lobes fade at the midline and rim. Explicit IDs exclude misplaced pelvic ligaments.
 # An illustration-guided soft-tissue contour, not a measured female muscle volume.
 gluteProjection=dict(partIds=['FJ1418','FJ1418M','FJ3513','FJ3606'],scope='BodyParts3D bilateral gluteus maximus and inferior gluteal veins only; posterior/inferior components excluded from bones, pelvic organs, ligaments and lumbar muscles',amplitude=.020,centerX=.070,centerY=.880,radiusX=.100,radiusY=.180,midlineRamp=[.030,.050],depthRamp=[-.120,-.095],lowerDepth=dict(heightRise=[.780,.860],heightFade=[.830,.895],depthRamp=[-.100,-.060]),inferior=dict(amplitude=.022,posteriorFade=[-.120,-.075],midlineRamp=[.006,.020],lateralFade=[.065,.105],heightRamp=[.760,.800,.840,.885],depthRamp=[-.090,-.060])),
 head=dict(scale=.95,center=[0.,1.61,-.03],ramp=[1.44,1.56]),
 # Lower nasal bridge and projection: depth beyond the face plane is compressed inside a smooth
 # window around the nose (nasal bones and cartilages), giving a flatter East Asian profile.
 nose=dict(plane=.062,scale=.52,center=[0.,1.578],radius=[.03,.045],rampDepth=.012),
)
def lateral_factor(y,knots=None):
    knots=MORPH['lateralKnots'] if knots is None else knots;out=np.full_like(y,knots[-1][1])
    out=np.where(y<knots[0][0],knots[0][1],out)
    for (y0,s0),(y1,s1) in zip(knots,knots[1:]):
        inside=(y>=y0)&(y<y1);out=np.where(inside,s0+(s1-s0)*smoothstep(y0,y1,y),out)
    return out
def morph(p,part_id=None):
    x,y,z=p[:,0],p[:,1],p[:,2];R=MORPH['lateralRadius']
    factor=lateral_factor(y);dx=(factor-1)*R*np.tanh(x/R)
    blend=MORPH.get('lateralBlend')
    if blend:
        reference=lateral_factor(y,blend['referenceKnots'])
        w=smoothstep(blend['innerRadius'],blend['outerRadius'],np.abs(x));r=blend['heightRamp']
        outer=np.sign(x)*blend['outerTranslation']*smoothstep(r[0],r[1],y)*(1-smoothstep(r[2],r[3],y))
        dx=(reference-1)*R*np.tanh(x/R)+(1-w)*(factor-reference)*R*np.tanh(x/R)+w*outer
    waist=MORPH.get('waistRefinement')
    if waist:
        r=waist['heightRamp'];weight=smoothstep(r[0],r[1],y)*(1-smoothstep(r[1],r[2],y))
        weight*=1-smoothstep(*waist['radialRamp'],np.abs(x))
        dx+=waist['delta']*R*np.tanh(x/R)*weight
    t=MORPH['thoraxDepth'];w=smoothstep(t['ramp'][0],t['ramp'][1],y)*(1-smoothstep(t['ramp'][2],t['ramp'][3],y))
    dz=z*(t['scale']-1)*w
    dy=np.zeros_like(y)
    g=MORPH.get('gluteProjection')
    if g and part_id in g['partIds']:
        r2=((np.abs(x)-g['centerX'])/g['radiusX'])**2+((y-g['centerY'])/g['radiusY'])**2
        posterior=1-smoothstep(*g['depthRamp'],z)
        lower=g.get('lowerDepth')
        if lower:
            low=(1-smoothstep(*lower['heightFade'],y))*(1-smoothstep(*lower['depthRamp'],z))
            if 'heightRise' in lower:low*=smoothstep(*lower['heightRise'],y)
            posterior=posterior+low-posterior*low
        gate=(1-smoothstep(0.,1.,r2))*smoothstep(*g['midlineRamp'],np.abs(x))*posterior
        dz-=g['amplitude']*gate
        inferior=g.get('inferior')
        if inferior:
            r=inferior['heightRamp']
            weight=smoothstep(*inferior['midlineRamp'],np.abs(x))*(1-smoothstep(*inferior['lateralFade'],np.abs(x)))
            weight*=smoothstep(r[0],r[1],y)*(1-smoothstep(r[2],r[3],y))*(1-smoothstep(*inferior['depthRamp'],z))
            if 'posteriorFade' in inferior:weight*=smoothstep(*inferior['posteriorFade'],z)
            dy-=inferior['amplitude']*weight
    q=np.stack([x+dx,y+dy,z+dz],axis=1)
    n=MORPH['nose'];r2=((x-n['center'][0])/n['radius'][0])**2+((y-n['center'][1])/n['radius'][1])**2
    wn=(1-smoothstep(0.,1.,r2))*smoothstep(n['plane']-n['rampDepth'],n['plane']+n['rampDepth'],z)
    q[:,2]=q[:,2]-wn*(1-n['scale'])*(z-n['plane'])
    h=MORPH['head'];wh=smoothstep(h['ramp'][0],h['ramp'][1],y)[:,None]
    q=q-wh*(1-h['scale'])*(q-np.array(h['center']))
    return q*MORPH['stature']
def morph_normals(p,n,part_id=None,h=1e-4):
    base=morph(p,part_id);J=np.empty((len(p),3,3))
    for axis in range(3):
        d=np.zeros(3);d[axis]=h;J[:,:,axis]=(morph(p+d,part_id)-base)/h
    cof=np.linalg.inv(J).transpose(0,2,1)
    out=np.einsum('nij,nj->ni',cof,n);return out/np.maximum(np.linalg.norm(out,axis=1,keepdims=True),1e-20)
# Jacobian sanity: sample the field over the body volume.
grid=np.stack(np.meshgrid(np.linspace(-.35,.35,29),np.linspace(0,1.75,71),np.linspace(-.15,.15,13),indexing='ij'),-1).reshape(-1,3)
# Screen both base field and the optional glute field. All eligible IDs use the same field.
min_dets={}
for label,part_id in [('shared',None),('glute','FJ1418')]:
    J=np.empty((len(grid),3,3));base=morph(grid,part_id)
    for axis in range(3):
        d=np.zeros(3);d[axis]=1e-4;J[:,:,axis]=(morph(grid+d,part_id)-base)/1e-4
    min_dets[label]=float(np.linalg.det(J).min());assert min_dets[label]>0,min_dets
min_det=min(min_dets.values())
for m in meshes:
    if m.get('group')=='mammary':continue
    part_id=m['part']['id']
    if part_id in MORPH['gluteProjection']['partIds']:
        m['provenance']['adaptation']+='; illustration-guided posterior and inferior glute contour (explicit part-ID scope)'
    m['nrm']=morph_normals(m['pos'],m['nrm'],part_id);m['pos']=morph(m['pos'],part_id)

# ---------------------------------------------------------------- 4. source-derived breast contour
breast_report=breast_contour.describe(fp)
breast_report['chestClearance']=dict(method='xy-grid z difference: adipose posterior (min z) minus pectoralis anterior (max z) on a 4 mm cell, inner 60% ellipse of each fat envelope. Positive is a front gap; negative is intersection. Geometric screening only.',pectoralis='muscular parts whose names match /pectoralis/i')
breast_report['sourceManifestSha256']=hashlib.sha256((INPUT/'atlas-female.json').read_bytes()).hexdigest()
breast_report['sourceChunks']={female['chunks'][i]['url']:hashlib.sha256(source_buffers[i]).hexdigest() for i in sorted({fp[key]['chunk'] for key in contoured})}
minimum_breast_jacobian=float('inf')
for m in meshes:
    if m.get('group')!='mammary':continue
    t=breast_report['transforms'][m['part']['id'][-1]]
    m['pos'],m['encodedNormals'],minimum=breast_contour.apply(m['pos'],m['nrm'],t)
    minimum_breast_jacobian=min(minimum_breast_jacobian,minimum)
breast_report['minimumVertexJacobian']=minimum_breast_jacobian

def _breast_chest_clearance(chunks_dir):
    return json.loads(subprocess.check_output(['node',str(ROOT/'scripts/breast-chest-clearance.mjs'),str(chunks_dir)],cwd=ROOT))

# ---------------------------------------------------------------- write atlas
for stale in glob.glob(str(OUT/'female-base-*.bin*')):os.remove(stale)
atlas=copy.deepcopy(male)
atlas.update(version='Female study prototype v3',sex='female',source='BodyParts3D framework + HRA female anatomy',scope='Estimated female proportions; experimental organ placement',reconstruction=True,parts=[],concepts=[],chunks=[])
chunks=atlas['chunks'];blob=bytearray()
def flush():
    global blob
    if not blob:return
    name=f'female-base-{len(chunks)}.bin';data=bytes(blob)
    compressed=gzip.compress(data,compresslevel=9,mtime=0)
    (OUT/name).write_bytes(data);(OUT/(name+'.gz')).write_bytes(compressed)
    chunks.append(dict(url=f'/models/{name}',bytes=len(data),gzip=f'/models/{name}.gz',gzipBytes=len(compressed)))
    blob=bytearray()
def append(a):
    while len(blob)%4:blob.append(0)
    offset=len(blob);blob.extend(a.tobytes());return offset
added=[]
for m in meshes:
    if len(blob)>4_000_000:flush()
    pos=m['pos'].astype('<f4');nrm=m['encodedNormals'] if 'encodedNormals' in m else np.rint(m['nrm']*32767).astype('<i2');idx=m['idx'].astype('<u4')
    p=copy.deepcopy(m['part'])
    if m.get('group')=='mammary' and any(k in p['id'] for k in ('nipple','areola')):p['system']='integumentary'
    p.update(chunk=len(chunks),vertexCount=len(pos),indexCount=len(idx),positions=append(pos),normals=append(nrm),indices=append(idx),bounds=[pos.min(axis=0).tolist(),pos.max(axis=0).tolist()],provenance=m['provenance'])
    atlas['parts'].append(p)
    if p['provenance']['source']!='BodyParts3D 4.0':p['group']=m['group'];added.append(p)
flush()

# Concepts must remain complete. Never relabel a truncated male body concept as female or
# leave hidden genital members in a compound selection.
retained_ids={p['id'] for p in retained}
for c in male['concepts']:
    if c['elements'] and all(e in retained_ids or e in replacements for e in c['elements']):
        atlas['concepts'].append(dict(c,elements=[r for e in c['elements'] for r in (replacements.get(e) or [e])]))
selected_ids={p['id'] for p in selected}
for c in female['concepts']:
    if c['elements'] and set(c['elements'])<=selected_ids:atlas['concepts'].append(copy.deepcopy(c))
for p in atlas['parts']:
    if not any(c['elements']==[p['id']] for c in atlas['concepts']):
        atlas['concepts'].append(dict(id='PART_'+p['id'],name=p['name'],elements=[p['id']]))
assert len({c['id'] for c in atlas['concepts']})==len(atlas['concepts'])
atlas['triangles']=sum(p['indexCount']//3 for p in atlas['parts'])
atlas['optimized']={'method':'Optimized BodyParts3D framework and fitted HRA non-breast anatomy follow the recorded body morph; original HRA breast topology follows a separate recorded contour in final coordinates','preservedMeshes':len(atlas['parts'])}
atlas.pop('sourceTriangles',None)
atlas['reconstructionReport']='/models/female-fit-report.json'

def extent(name,axis):
    ps=[p for p in atlas['parts'] if p['name']==name];src=[mp[p['id']] for p in ps]
    return dict(before=round(max(q['bounds'][1][axis] for q in src)-min(q['bounds'][0][axis] for q in src),4),after=round(max(q['bounds'][1][axis] for q in ps)-min(q['bounds'][0][axis] for q in ps),4))
def span(names,axis,after_names=None):
    src=[p for p in male['parts'] if p['name'] in names];ps=[p for p in atlas['parts'] if p['name'] in (after_names or names)]
    return dict(before=round(max(q['bounds'][1][axis] for q in src)-min(q['bounds'][0][axis] for q in src),4),after=round(max(q['bounds'][1][axis] for q in ps)-min(q['bounds'][0][axis] for q in ps),4))
landmarks=dict(
 stature=dict(before=round(max(q['bounds'][1][1] for q in mp.values()),4),after=round(max(p['bounds'][1][1] for p in atlas['parts']),4)),
 biacromialWidth=span(['Left scapula','Right scapula'],0),
 biIliacWidth=span(['Left hip bone','Right hip bone'],0,['Left ilium','Right ilium']),
 headWidth=span(['Left parietal bone','Right parietal bone'],0),
 chestDepth=extent('Body of sternum',2),
)
added_entries=[dict(id=p['id'],name=p['name'],system=p['system'],transform='breastContour' if p['group']=='mammary' else p['group']) for p in added]
for p in atlas['parts']:p.pop('group',None)
(OUT/'atlas-female-reconstructed.json').write_text(json.dumps(atlas,separators=(',',':')))
clearance=_breast_chest_clearance(OUT)
report=dict(method='BodyParts3D framework and fitted HRA non-breast structures follow the body morph; 16 original HRA breast meshes receive a separate contour and final-frame translation',reviewStatus='Experimental; proportions estimated, organ placement unreviewed',transforms=transforms,morph=MORPH,breastContour=breast_report,regenerated=[],contoured=contoured,replacements=replacements,landmarks=landmarks,retained=[p['id'] for p in retained],added=added_entries,excluded=excluded,checks=dict(retainedGeometryUnchanged=False,limbProportionsUnchanged=False,minimumMorphJacobian=min_det,minimumMorphJacobianByField=min_dets,minimumTransformDeterminant=float(min(np.prod(t['scale']) for t in transforms.values())),minimumBreastVertexJacobian=minimum_breast_jacobian,breastSourceTopologyPreserved=True,breastChestCoreCoveredCells=clearance['coreCoveredCells'],breastChestCoreGapMinM=clearance['coreGapMinM'],breastChestCoreGapMedianM=clearance['coreGapMedianM'],breastChestCoreGapMaxM=clearance['coreGapMaxM']),parts=len(atlas['parts']),concepts=len(atlas['concepts']),triangles=atlas['triangles'])
(OUT/'female-fit-report.json').write_text(json.dumps(report,indent=2))
# The published ligament recipe is evaluated against the freshly rebuilt base.
# The helper verifies source pins, the joint field and unchanged base buffers.
with tempfile.TemporaryDirectory(prefix='female-joint-rebuild-') as temporary:
    staged=Path(temporary)
    subprocess.run([sys.executable,str(ROOT/'scripts/publish-female-joint-additions.py'),
                    '--candidate-dir',str(ROOT/'data/anatomy/joint-additions'),
                    '--base-models-dir',str(OUT),'--source-models-dir',str(INPUT),
                    '--output-dir',str(staged)],check=True)
    enriched=json.loads((staged/'atlas-female-reconstructed.json').read_text())
    for chunk in enriched['chunks'][len(chunks):]:
        for field in ['url','gzip']:
            name=Path(chunk[field]).name;shutil.copyfile(staged/name,OUT/name)
    for name in ['atlas-female-reconstructed.json','female-fit-report.json']:
        shutil.copyfile(staged/name,OUT/name)
    atlas=enriched;chunks=atlas['chunks'];report=json.loads((OUT/'female-fit-report.json').read_text())
print(json.dumps({k:report[k] for k in ['checks','landmarks','parts','concepts','triangles']},indent=2))
print('Retained',len(retained),'base meshes; added',len(added)+len(report['jointAdditions']['ids']),'female meshes; excluded',len(excluded),'base meshes')
print('Compressed MB',sum(c['gzipBytes'] for c in chunks)/1e6)
