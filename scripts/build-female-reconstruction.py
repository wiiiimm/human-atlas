"""Build a female study prototype from the BodyParts3D framework and HRA female organs.

Run: python3 scripts/build-female-reconstruction.py (requires numpy).

Pipeline
1. Omit male-specific BodyParts3D structures from the manifest.
2. Fit HRA female reproductive organs and breast tissue into the base pelvis and chest.
3. Drape the breast assembly onto the retained chest wall so it rests on the pectoral muscles.
4. Apply one smooth whole-body morph (stature, shoulders, thorax, waist, pelvis, head) to
   every mesh, so bones, muscles, vessels, and fitted organs deform together and stay attached.

The morph is a continuous displacement field with a positive Jacobian everywhere, so no mesh
tears and attachments are preserved. Its parameters are recorded in the fit report and are
re-applied by scripts/validate-atlas.mjs to verify every vertex. This is a study prototype
with estimated female proportions, not an anatomically validated female atlas.
"""
from pathlib import Path
import copy, glob, gzip, json, os, re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'public/models'
male=json.loads((OUT/'atlas.json').read_text())
female=json.loads((OUT/'atlas-female.json').read_text())
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
# Breast assembly: keep HRA size, place the nipple line near the fourth intercostal space of the
# base chest (about 71% of stature, matching its relative height in the HRA body).
breast_scale=np.array([1.03,1.,1.])
breast_offset=np.array([.009,.045,.115])
# Female pelvis: hip bone envelope matched to the male hip bone envelope so the acetabula meet
# the retained femoral heads and the iliac crests keep their height.
hip_ids=[i for k in ('FJ3152','FJ3288') for i in replacements[k]]
pelvis_scale,pelvis_offset=bounds_fit([fp[i] for i in hip_ids],[mp['FJ3152'],mp['FJ3288']])
transforms={
 'reproductive':dict(scale=pelvic_scale.tolist(),offset=pelvic_offset.tolist(),basis='HRA bladder bounds aligned to retained BodyParts3D bladder bounds'),
 'mammary':dict(scale=breast_scale.tolist(),offset=breast_offset.tolist(),basis='HRA breast size retained; nipple line placed at the base chest fourth intercostal level, then draped onto the chest wall'),
 'pelvis':dict(scale=pelvis_scale.tolist(),offset=pelvis_offset.tolist(),basis='HRA hip bone bounds aligned to BodyParts3D hip bone bounds; sacrum and coccyx share the transform')
}
pelvis_ids=[i for ids in replacements.values() for i in ids]
selected=[p for p in female['parts'] if p['system']=='reproductive' or (p['system']=='integumentary' and p['id']!='VH_F_skin') or p['id'] in pelvis_ids]

# Assemble every mesh in the base coordinate frame before deformation.
meshes=[]  # dict(part, pos, nrm, idx, provenance)
for p in retained:
    pos,nrm,idx=read(male_buffers,p)
    meshes.append(dict(part=p,pos=pos,nrm=nrm,idx=idx,provenance={'source':'BodyParts3D 4.0','sourceId':p['id'],'adaptation':'Retained reference geometry; reshaped by the shared female body morph'}))
for source in selected:
    group='pelvis' if source['id'] in pelvis_ids else 'mammary' if source['system']=='integumentary' else 'reproductive'
    t=transforms[group];scale=np.array(t['scale']);offset=np.array(t['offset'])
    pos,nrm,idx=read(source_buffers,source)
    pos=pos*scale+offset
    nrm=nrm/scale;nrm/=np.maximum(np.linalg.norm(nrm,axis=1,keepdims=True),1e-20)
    p=copy.deepcopy(source);p['system']='skeletal' if group=='pelvis' else group
    if group=='pelvis':p['name']=pelvis_names[source['id']]
    meshes.append(dict(part=p,pos=pos,nrm=nrm,idx=idx,group=group,provenance={'source':'HRA united-female v1.5','sourceId':source['id'],'adaptation':('Female pelvis fitted to the BodyParts3D hip bone envelope, replacing the male pelvis' if group=='pelvis' else 'Female reference anatomy fitted to the BodyParts3D framework'+(' and draped onto the chest wall' if group=='mammary' else ''))+'; experimental placement reshaped by the shared female body morph'}))

# ---------------------------------------------------------------- 3. breast drape
# Build a depth map of the anterior chest wall (front-most z of chest muscles and bones on a
# regular x/y grid), a matching map of the breast assembly's posterior envelope, and shift each
# breast column in z so the assembly rests a few millimetres into the wall instead of floating.
CELL=.004;GX0,GX1,GY0,GY1=-.24,.24,1.02,1.50
NX=int(round((GX1-GX0)/CELL))+1;NY=int(round((GY1-GY0)/CELL))+1
def cell_index(pos):
    ix=np.clip(np.rint((pos[:,0]-GX0)/CELL).astype(int),0,NX-1);iy=np.clip(np.rint((pos[:,1]-GY0)/CELL).astype(int),0,NY-1);return ix,iy
def project(items,reducer):
    # Rasterize every triangle onto the grid (barycentric depth at cell centres), so cells between
    # the sparse vertices of simplified meshes are covered too. reducer is max (front) or min (back).
    grid=np.full((NY,NX),np.nan);pick=np.fmax if reducer is max else np.fmin
    for m in items:
        pos=m['pos'];tri=pos[m['idx'].reshape(-1,3)]
        for a,b,c in tri:
            if max(a[1],b[1],c[1])<GY0 or min(a[1],b[1],c[1])>GY1 or max(a[0],b[0],c[0])<GX0 or min(a[0],b[0],c[0])>GX1:continue
            i0=max(0,int(np.ceil((min(a[0],b[0],c[0])-GX0)/CELL)));i1=min(NX-1,int(np.floor((max(a[0],b[0],c[0])-GX0)/CELL)))
            j0=max(0,int(np.ceil((min(a[1],b[1],c[1])-GY0)/CELL)));j1=min(NY-1,int(np.floor((max(a[1],b[1],c[1])-GY0)/CELL)))
            if i1<i0 or j1<j0:continue
            xs=GX0+np.arange(i0,i1+1)*CELL;ys=GY0+np.arange(j0,j1+1)*CELL;X,Y=np.meshgrid(xs,ys)
            det=(b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1])
            if abs(det)<1e-12:continue
            u=((X-a[0])*(c[1]-a[1])-(c[0]-a[0])*(Y-a[1]))/det;v=((b[0]-a[0])*(Y-a[1])-(X-a[0])*(b[1]-a[1]))/det
            inside=(u>=-1e-9)&(v>=-1e-9)&(u+v<=1+1e-9)
            if not inside.any():continue
            z=a[2]+u*(b[2]-a[2])+v*(c[2]-a[2])
            sub=grid[j0:j1+1,i0:i1+1];sub[inside]=pick(sub[inside],z[inside])
    # Vertices still count, so tiny meshes thinner than a cell are not lost.
    for m in items:
        pos=m['pos'];ix,iy=cell_index(pos)
        keep=(pos[:,0]>=GX0)&(pos[:,0]<=GX1)&(pos[:,1]>=GY0)&(pos[:,1]<=GY1)
        np.fmax.at(grid,(iy[keep],ix[keep]),pos[keep,2]) if reducer is max else np.fmin.at(grid,(iy[keep],ix[keep]),pos[keep,2])
    return grid
def fill_nearest(grid):
    filled=grid.copy();mask=~np.isnan(filled)
    while not mask.all():
        grown=filled.copy()
        for dy,dx in [(1,0),(-1,0),(0,1),(0,-1)]:
            shifted=np.roll(filled,(dy,dx),axis=(0,1));shifted_mask=np.roll(mask,(dy,dx),axis=(0,1))
            take=(~mask)&shifted_mask&np.isnan(grown);grown[take]=shifted[take]
        filled=grown;mask=~np.isnan(filled)
    return filled
def blur(grid,sigma):
    radius=int(3*sigma);k=np.exp(-np.arange(-radius,radius+1)**2/(2*sigma**2));k/=k.sum()
    out=np.apply_along_axis(lambda r:np.convolve(np.pad(r,radius,mode='edge'),k,'valid'),1,grid)
    return np.apply_along_axis(lambda c:np.convolve(np.pad(c,radius,mode='edge'),k,'valid'),0,out)
chest=[m for m in meshes if m['part']['system'] in('muscular','skeletal','connective') and m['part']['bounds'][1][2]>.02 and m['part']['bounds'][0][1]<GY1 and m['part']['bounds'][1][1]>GY0 and abs(m['part']['bounds'][0][0]+m['part']['bounds'][1][0])/2<.26 and 'papillary' not in m['part']['name'] and m['part']['name']!='Diaphragm']
wall=blur(fill_nearest(project(chest,max)),1.)
breast=[m for m in meshes if m['part']['system']=='mammary']
envelope=project(breast,min);covered=~np.isnan(envelope)
# A small standoff keeps the front of the thin upper breast clear of the muscle after smoothing;
# the max filter lets "push forward" win locally so the pectoral ridge never pokes through.
EMBED=-.002
shift=np.where(covered,wall-EMBED-envelope,np.nan)
def max_filter(grid,radius):
    out=grid.copy()
    for dy in range(-radius,radius+1):
        for dx in range(-radius,radius+1):
            out=np.fmax(out,np.roll(grid,(dy,dx),axis=(0,1)))
    return out
drape_grid=blur(max_filter(fill_nearest(shift),2),1.5)
def bilinear(grid,pos):
    fx=np.clip((pos[:,0]-GX0)/CELL,0,NX-1.000001);fy=np.clip((pos[:,1]-GY0)/CELL,0,NY-1.000001)
    x0=np.floor(fx).astype(int);y0=np.floor(fy).astype(int);tx=fx-x0;ty=fy-y0
    return (grid[y0,x0]*(1-tx)*(1-ty)+grid[y0,x0+1]*tx*(1-ty)+grid[y0+1,x0]*(1-tx)*ty+grid[y0+1,x0+1]*tx*ty)
for m in breast:
    m['pos']=m['pos'].copy();m['pos'][:,2]+=bilinear(drape_grid,m['pos'])
drape=dict(cell=CELL,origin=[GX0,GY0],columns=NX,rows=NY,embed=EMBED,values=np.round(drape_grid,6).tolist())

# Regenerate each breast fat body as a heightfield that grows out of the chest wall. The HRA
# interlobar fat keeps its real thickness profile (front minus back per column), but its edge is
# feathered over ~2.5 cm and dips below the wall, so the mound rises from the pectoral surface
# instead of sitting on it as a separate blob. The glandular structures stay inside it.
FEATHER=6;STANDOFF=.003;UNDERCUT=.002
def closing(mask,r):
    grow=mask.copy()
    for _ in range(r):grow=grow|np.roll(grow,1,0)|np.roll(grow,-1,0)|np.roll(grow,1,1)|np.roll(grow,-1,1)
    for _ in range(r):grow=grow&np.roll(grow,1,0)&np.roll(grow,-1,0)&np.roll(grow,1,1)&np.roll(grow,-1,1)
    return grow
def distance_outside(mask,limit):
    dist=np.where(mask,0,np.inf);frontier=mask.copy()
    for step in range(1,limit+1):
        grown=frontier|np.roll(frontier,1,0)|np.roll(frontier,-1,0)|np.roll(frontier,1,1)|np.roll(frontier,-1,1)
        dist=np.where(grown&~frontier,step,dist);frontier=grown
    return dist
def heightfield_mesh(mask,front,back):
    full=mask[:-1,:-1]&mask[1:,:-1]&mask[:-1,1:]&mask[1:,1:]
    xs=GX0+np.arange(NX)*CELL;ys=GY0+np.arange(NY)*CELL
    fid=np.full((NY,NX),-1);bid=np.full((NY,NX),-1);verts=[]
    def vid(table,j,i,z):
        if table[j,i]<0:table[j,i]=len(verts);verts.append([xs[i],ys[j],z[j,i]])
        return table[j,i]
    tris=[]
    for j,i in zip(*np.nonzero(full)):
        a,b,c,d=[(j,i),(j,i+1),(j+1,i+1),(j+1,i)]
        fa,fb,fc,fd=[vid(fid,*q,front) for q in (a,b,c,d)];ba,bb,bc,bd=[vid(bid,*q,back) for q in (a,b,c,d)]
        tris+=[[fa,fb,fc],[fa,fc,fd],[ba,bc,bb],[ba,bd,bc]]
    padded=np.pad(full,1)
    for j in range(NY):
        for i in range(NX-1):  # edge (j,i)-(j,i+1): cells below (j-1) and above (j)
            above=padded[j+1,i+1];below=padded[j,i+1]
            if above==below:continue
            f0,f1=vid(fid,j,i,front),vid(fid,j,i+1,front);b0,b1=vid(bid,j,i,back),vid(bid,j,i+1,back)
            quad=[[f0,f1,b1],[f0,b1,b0]]  # normal +y
            tris+=quad if below and not above else [t[::-1] for t in quad]
    for i in range(NX):
        for j in range(NY-1):  # edge (j,i)-(j+1,i): cells left (i-1) and right (i)
            right=padded[j+1,i+1];left=padded[j+1,i]
            if right==left:continue
            f0,f1=vid(fid,j,i,front),vid(fid,j+1,i,front);b0,b1=vid(bid,j,i,back),vid(bid,j+1,i,back)
            quad=[[f0,b1,f1],[f0,b0,b1]]  # normal +x
            tris+=quad if left and not right else [t[::-1] for t in quad]
    pos=np.array(verts);idx=np.array(tris,dtype=np.uint32)
    n=np.zeros_like(pos);fn=np.cross(pos[idx[:,1]]-pos[idx[:,0]],pos[idx[:,2]]-pos[idx[:,0]])
    for k in range(3):np.add.at(n,idx[:,k],fn)
    n/=np.maximum(np.linalg.norm(n,axis=1,keepdims=True),1e-20)
    return pos,n,idx.reshape(-1)
regenerated=[]
for m in breast:
    if 'Interlobar adipose' not in m['part']['name']:continue
    side='left' if m['part']['id'].endswith('_L') else 'right'
    # Envelope of the whole draped breast assembly on this side (the interlobar fat alone has
    # cavities where the lobes sit). The nipple and areolar tubercles stay proud of the surface.
    assembly=[dict(pos=b['pos'],idx=b['idx']) for b in breast if b['part']['id'].endswith('_L' if side=='left' else '_R')]
    for b,src in zip(assembly,[b for b in breast if b['part']['id'].endswith('_L' if side=='left' else '_R')]):b['part']=src['part']
    outer=[b for b in assembly if not any(k in b['part']['id'] for k in ('nipple','areolar_tubercles'))]
    F=project(outer,max);B=project(assembly,min);footprint=closing(~np.isnan(F),2)
    T=np.where(~np.isnan(F),F-B,np.nan);T=np.where(footprint,fill_nearest(T),0.)
    T=blur(max_filter(np.where(footprint,T,0.),2),1.5)+.002
    d=distance_outside(footprint,FEATHER);feather=1-smoothstep(0,FEATHER,np.minimum(d,FEATHER))
    T=T*feather;standoff=STANDOFF*feather-UNDERCUT*(1-feather)
    front=wall+T+standoff;back=wall-UNDERCUT
    mask=(front-wall>=-UNDERCUT*.5)&(d<=FEATHER)
    pos,nrm,idx=heightfield_mesh(mask,front,back)
    m.update(pos=pos,nrm=nrm,idx=idx)
    m['part']['name']=f'Adipose tissue of {side} breast'
    m['provenance']['adaptation']='Regenerated from the HRA interlobar adipose thickness profile as a feathered heightfield on the chest wall; reshaped by the shared female body morph'
    regenerated.append(m['part']['id'])
gap_before=float(np.nanmedian(shift));residual=np.where(covered,wall-EMBED-(envelope+drape_grid),np.nan);gap_after=float(np.nanmedian(residual));gap_max=float(np.nanmax(residual))

# ---------------------------------------------------------------- 4. whole-body morph
# Lateral width factor by height (base frame, before stature scaling). Between knots the
# factor follows a smoothstep so the field is continuous and slopes stay gentle along the arms.
MORPH=dict(
 stature=.95,
 lateralKnots=[[0.,1.],[.55,1.],[.80,1.05],[.92,1.10],[1.03,1.10],[1.21,.96],[1.29,.95],[1.42,.91],[1.48,.91],[1.60,1.]],
 lateralRadius=.17,  # beyond this |x| the lateral change saturates into a translation, so arms move with the torso instead of being squeezed
 thoraxDepth=dict(scale=.96,ramp=[1.05,1.15,1.35,1.45]),
 head=dict(scale=.95,center=[0.,1.61,-.03],ramp=[1.44,1.56]),
)
def lateral_factor(y):
    knots=MORPH['lateralKnots'];out=np.full_like(y,knots[-1][1])
    out=np.where(y<knots[0][0],knots[0][1],out)
    for (y0,s0),(y1,s1) in zip(knots,knots[1:]):
        inside=(y>=y0)&(y<y1);out=np.where(inside,s0+(s1-s0)*smoothstep(y0,y1,y),out)
    return out
def morph(p):
    x,y,z=p[:,0],p[:,1],p[:,2];R=MORPH['lateralRadius']
    dx=np.sign(x)*(lateral_factor(y)-1)*R*np.tanh(np.abs(x)/R)
    t=MORPH['thoraxDepth'];w=smoothstep(t['ramp'][0],t['ramp'][1],y)*(1-smoothstep(t['ramp'][2],t['ramp'][3],y))
    dz=z*(t['scale']-1)*w
    q=np.stack([x+dx,y,z+dz],axis=1)
    h=MORPH['head'];wh=smoothstep(h['ramp'][0],h['ramp'][1],y)[:,None]
    q=q-wh*(1-h['scale'])*(q-np.array(h['center']))
    return q*MORPH['stature']
def morph_normals(p,n,h=1e-4):
    base=morph(p);J=np.empty((len(p),3,3))
    for axis in range(3):
        d=np.zeros(3);d[axis]=h;J[:,:,axis]=(morph(p+d)-base)/h
    cof=np.linalg.inv(J).transpose(0,2,1)
    out=np.einsum('nij,nj->ni',cof,n);return out/np.maximum(np.linalg.norm(out,axis=1,keepdims=True),1e-20)
# Jacobian sanity: sample the field over the body volume.
grid=np.stack(np.meshgrid(np.linspace(-.35,.35,29),np.linspace(0,1.75,71),np.linspace(-.15,.15,13),indexing='ij'),-1).reshape(-1,3)
J=np.empty((len(grid),3,3));base=morph(grid)
for axis in range(3):
    d=np.zeros(3);d[axis]=1e-4;J[:,:,axis]=(morph(grid+d)-base)/1e-4
min_det=float(np.linalg.det(J).min());assert min_det>0,min_det

for m in meshes:
    m['nrm']=morph_normals(m['pos'],m['nrm']);m['pos']=morph(m['pos'])

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
    pos=m['pos'].astype('<f4');nrm=np.rint(m['nrm']*32767).astype('<i2');idx=m['idx'].astype('<u4')
    p=copy.deepcopy(m['part'])
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
atlas['optimized']={'method':'Optimized BodyParts3D framework and affine-fitted optimized HRA female meshes, reshaped by one smooth female body morph','preservedMeshes':len(atlas['parts'])}
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
report=dict(method='BodyParts3D framework with fitted and draped HRA female structures, reshaped by one smooth whole-body morph',reviewStatus='Experimental; proportions estimated, organ placement unreviewed',transforms=transforms,morph=MORPH,drape=drape,regenerated=regenerated,replacements=replacements,landmarks=landmarks,retained=[p['id'] for p in retained],added=[dict(id=p['id'],name=p['name'],system=p['system'],transform=p['group']) for p in added],excluded=excluded,checks=dict(retainedGeometryUnchanged=False,limbProportionsUnchanged=False,minimumMorphJacobian=min_det,minimumTransformDeterminant=float(min(np.prod(t['scale']) for t in transforms.values())),breastWallGapBeforeDrapeM=gap_before,breastWallGapAfterDrapeM=gap_after,breastWallMaxResidualM=gap_max),parts=len(atlas['parts']),concepts=len(atlas['concepts']),triangles=atlas['triangles'])
for p in atlas['parts']:p.pop('group',None)
(OUT/'atlas-female-reconstructed.json').write_text(json.dumps(atlas,separators=(',',':')))
(OUT/'female-fit-report.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:report[k] for k in ['checks','landmarks','parts','concepts','triangles']},indent=2))
print('Retained',len(retained),'base meshes; added',len(added),'female meshes; excluded',len(excluded),'base meshes')
print('Compressed MB',sum(c['gzipBytes'] for c in chunks)/1e6)
