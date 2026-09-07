#!/usr/bin/env python3
"""Small source/candidate fixtures test append-only packaging and rejection gates."""
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

spec = importlib.util.spec_from_file_location('append_joints', Path(__file__).with_name('publish-female-joint-additions.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class JointAdditionsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='joint-additions-test-')
        self.folder = Path(self.temp.name); self.root = self.folder/'repo'
        self.models = self.root/'public/models'; self.models.mkdir(parents=True)
        self.candidate = self.folder/'candidate'; self.candidate.mkdir()
        self.output = self.folder/'output'
        vertices = np.array([[.10,.40,.02],[.11,.40,.02],[.10,.41,.02],[.10,.40,.03]], dtype='<f4')
        faces = np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]], dtype='<u4')
        self.positions = vertices; self.faces = faces
        def model(ids, name):
            payload = bytearray(); parts = []
            for key in ids:
                offsets = []
                for array in [vertices, np.rint(m.area_normals(vertices, faces)*32767).astype('<i2'), faces]:
                    payload.extend(b'\0'*((-len(payload)) % 4)); offsets.append(len(payload)); payload.extend(array.tobytes())
                parts.append(dict(id=key, name='anterolateral ligament of knee' if key in m.IDS else 'Existing part', system='skeletal', conceptId='SOURCE:ONE', chunk=0, positions=offsets[0], normals=offsets[1], indices=offsets[2], vertexCount=4, indexCount=12, bounds=[vertices.min(0).tolist(), vertices.max(0).tolist()]))
            data=bytes(payload); compressed=gzip.compress(data,mtime=0)
            (self.models/(name+'.bin')).write_bytes(data); (self.models/(name+'.bin.gz')).write_bytes(compressed)
            return dict(parts=parts, concepts=[dict(id='BASE_CONCEPT',name='Existing concept',elements=[ids[0]])], chunks=[dict(url='/models/'+name+'.bin',bytes=len(data),gzip='/models/'+name+'.bin.gz',gzipBytes=len(compressed))], triangles=4*len(ids), optimized=dict(preservedMeshes=len(ids)), custom='Preserve this metadata')
        self.base = model(['BASE'], 'base'); source=model(m.IDS, 'source')
        morph=dict(lateralKnots=[[0,1],[2,1]],lateralRadius=.1,thoraxDepth=dict(scale=1,ramp=[0,.5,1.5,2]),nose=dict(center=[0,2],radius=[1,1],plane=0,rampDepth=1,scale=1),head=dict(scale=1,ramp=[0,2],center=[0,0,0]),stature=1)
        self.fit=dict(morph=morph,parts=1,concepts=1,triangles=4,added=[{'id':'EXISTING_IMPORT'}],retained=['BASE'],checks={'original':'unchanged'})
        for name,data in [('atlas-female-reconstructed.json',self.base),('atlas-female.json',source),('female-fit-report.json',self.fit)]:
            (self.models/name).write_bytes(m.encoded(data))
        hashes={str(p.relative_to(self.root)):m.sha(p.read_bytes()) for p in self.models.iterdir()}
        report=dict(existingAtlasByteIdentityVerified=True,inputHashes=hashes,method={'registration':'Fixture identity'},sides={})
        npz_hashes={}
        for key in m.IDS:
            normals=m.area_normals(vertices,faces).astype('<f4')
            path=self.candidate/(key+'.npz');np.savez(path,positions=vertices,normals=normals,indices=faces)
            npz_hashes[key]=m.sha(path.read_bytes())
            report['sides'][key[-1]]=dict(fit=dict(scale=1,rotationRowVectors=np.eye(3).tolist(),translationM=[0,0,0]),refinement=dict(ridge=.08,kernelRadiusM=.025,controlSourcePointsM=[[0,0,0]],coefficients=[[0,0,0]]),parts=[dict(id=key,topologyPreserved=True,distortion=dict(nonpositiveJacobianSamples=0,minimumJacobianDeterminant=1),export=dict(npz=path.name,positionsSha256=m.sha(vertices.tobytes()),indicesSha256=m.sha(faces.tobytes()),vertices=4,triangles=4,boundsM=[vertices.min(0).tolist(),vertices.max(0).tolist()]))])
        (self.candidate/'report.json').write_bytes(m.encoded(report))
        self.report=report
        self.patches=[patch.object(m,'APPROVED_REPORT_SHA',m.sha((self.candidate/'report.json').read_bytes())),patch.object(m,'APPROVED_NPZ_SHA',npz_hashes)]
        for item in self.patches:item.start()

    def tearDown(self):
        for item in reversed(self.patches):item.stop()
        self.temp.cleanup()

    def build(self, output=None, base=None):
        return m.run(self.candidate,output or self.output,root=self.root,base_models_dir=base,source_models_dir=self.models)

    def repin_report(self):
        (self.candidate/'report.json').write_bytes(m.encoded(self.report))
        m.APPROVED_REPORT_SHA=m.sha((self.candidate/'report.json').read_bytes())

    def replace_array(self,key,array_name,value):
        path=self.candidate/(key+'.npz')
        with np.load(path) as archive:data={name:archive[name] for name in archive.files}
        data[array_name]=value;np.savez(path,**data)
        m.APPROVED_NPZ_SHA[key]=m.sha(path.read_bytes())
        if array_name in ['positions','indices']:
            self.report['sides'][key[-1]]['parts'][0]['export'][array_name+'Sha256']=m.sha(value.tobytes())
        self.repin_report()

    def test_exact_append_only_deterministic_payload_and_metadata(self):
        manifest,fit=self.build()
        self.assertEqual(manifest['parts'][:-2],self.base['parts'])
        self.assertEqual(manifest['concepts'][:-2],self.base['concepts'])
        self.assertEqual(manifest['chunks'][:-1],self.base['chunks'])
        base_on_disk=json.loads((self.models/'atlas-female-reconstructed.json').read_text())
        output_on_disk=json.loads((self.output/'atlas-female-reconstructed.json').read_text())
        for key in ['parts','concepts','chunks']:
            count=len(base_on_disk[key])
            self.assertEqual(json.dumps(output_on_disk[key][:count],separators=(',',':')),json.dumps(base_on_disk[key],separators=(',',':')))
        self.assertEqual({k:v for k,v in manifest.items() if k not in ['parts','concepts','chunks','triangles','jointAdditions','optimized']},{k:v for k,v in self.base.items() if k not in ['parts','concepts','chunks','triangles','optimized']})
        self.assertEqual({k:v for k,v in fit.items() if k not in ['parts','concepts','triangles','jointAdditions']},{k:v for k,v in self.fit.items() if k not in ['parts','concepts','triangles']})
        self.assertEqual(manifest['optimized'], dict(self.base['optimized'],preservedMeshes=len(manifest['parts'])))
        self.assertEqual([p['id'] for p in manifest['parts'][-2:]],list(m.IDS))
        for p in manifest['parts'][-2:]:
            positions,normals,faces=m.arrays(self.output,manifest,p)
            np.testing.assert_array_equal(positions,self.positions);np.testing.assert_array_equal(faces,self.faces)
            self.assertEqual(p['system'],'connective');self.assertTrue(np.all(abs(np.linalg.norm(normals.astype(float)/32767,axis=1)-1)<5e-5))
            self.assertEqual(p['positions']%4,0);self.assertEqual(p['indices']%4,0)
        before={p.name:m.sha(p.read_bytes()) for p in self.output.iterdir()}
        self.build();self.assertEqual(before,{p.name:m.sha(p.read_bytes()) for p in self.output.iterdir()})
        for chunk in manifest['chunks']:
            self.assertEqual(gzip.decompress((self.output/Path(chunk['gzip']).name).read_bytes()),(self.output/Path(chunk['url']).name).read_bytes())
        self.assertEqual((self.output/'base.bin').read_bytes(),(self.models/'base.bin').read_bytes())

    def test_stale_report_and_changed_source_are_rejected_before_writes(self):
        path=self.candidate/'report.json';path.write_bytes(path.read_bytes()+b' ')
        with self.assertRaisesRegex(ValueError,'pinned reviewed'):self.build()
        self.assertFalse(self.output.exists())
        self.repin_report();(self.models/'source.bin').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'Stale candidate input'):self.build()
        self.assertFalse(self.output.exists())

    def test_rehashed_topology_mutation_cannot_bypass_source_gate(self):
        self.replace_array(m.IDS[0],'indices',self.faces[:,::-1].copy())
        with self.assertRaisesRegex(ValueError,'triangle indices changed'):self.build()
        self.assertFalse(self.output.exists())

    def test_rehashed_position_mutation_cannot_bypass_recorded_transform(self):
        self.replace_array(m.IDS[0],'positions',self.positions+np.float32(.001))
        with self.assertRaisesRegex(ValueError,'recorded source transform'):self.build()

    def test_nonunit_normals_and_unselected_candidate_config_are_rejected(self):
        self.replace_array(m.IDS[0],'normals',np.zeros_like(self.positions))
        with self.assertRaisesRegex(ValueError,'not unit'):self.build()
        self.replace_array(m.IDS[0],'normals',m.area_normals(self.positions,self.faces).astype('<f4'))
        self.report['sides']['L']['refinement']['ridge']=.02;self.repin_report()
        with self.assertRaisesRegex(ValueError,'Only reviewed'):self.build()

    def test_external_staged_base_and_output_symlink_guards(self):
        staged=self.folder/'staged';staged.mkdir()
        for path in self.models.iterdir(): (staged/path.name).write_bytes(path.read_bytes())
        self.build(output=staged,base=staged)
        self.assertEqual(json.loads((staged/'atlas-female-reconstructed.json').read_text())['parts'][:1],self.base['parts'])
        linked=self.folder/'linked';linked.symlink_to(self.models,target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'outside'):self.build(output=linked)
        self.output.mkdir();(self.output/'atlas-female-reconstructed.json').symlink_to(self.models/'atlas-female-reconstructed.json')
        with self.assertRaisesRegex(ValueError,'symlinks'):self.build()
        self.assertEqual(json.loads((self.models/'atlas-female-reconstructed.json').read_text()),self.base)


if __name__=='__main__':unittest.main()
