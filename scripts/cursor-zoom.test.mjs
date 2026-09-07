import fs from 'node:fs';
import assert from 'node:assert/strict';
import * as T from 'three';
import ts from 'typescript';
import {partIsVisible} from '../app/anatomy.ts';
const source=fs.readFileSync(new URL('../app/scene.tsx',import.meta.url),'utf8');
const start=source.indexOf('  const wheelPlane='),end=source.indexOf("  renderer.domElement.addEventListener('wheel'",start);
assert(start>0&&end>start);
const compiled=ts.transpileModule(source.slice(start,end)+'\nreturn wheel;', {compilerOptions:{target:ts.ScriptTarget.ES2022,module:ts.ModuleKind.None}}).outputText;
globalThis.WheelEvent={DOM_DELTA_PIXEL:0,DOM_DELTA_LINE:1,DOM_DELTA_PAGE:2};
const rect={left:30,top:20,width:1000,height:800};
const make=(specs=[],options={})=>{
 const camera=new T.PerspectiveCamera(34,rect.width/rect.height,.005,100);camera.position.set(0,0,4);camera.lookAt(0,0,0);if(options.viewOffset)camera.setViewOffset(1000,800,130,-80,1000,800);camera.updateMatrixWorld();
 const controls={enabled:true,enableZoom:true,zoomSpeed:1,minDistance:.07,maxDistance:40,target:new T.Vector3(),...options.controls};
 const atlas={parts:specs.map((s,i)=>({id:s.id??String(i),system:s.system??'muscular'}))};
 const bounds=specs.map(s=>new T.Box3(new T.Vector3(-1,-1,s.z-.1),new T.Vector3(1,1,s.z+.1)));
 const pickers=specs.map((s,i)=>{const g=new T.BoxGeometry(2,2,.2);g.translate(0,0,s.z);const mesh=new T.Mesh(g);mesh.updateMatrixWorld();return mesh;});
 const data=new Float32Array(specs.length*4);specs.forEach((s,i)=>data[i*4+3]=s.rendered===false?0:1);
 const latest={current:{visible:['muscular','integumentary','mammary'],selected:[],isolate:false,breastView:'tissue',...options.state}};
 const hover={hidden:false},renderer={domElement:{getBoundingClientRect:()=>rect}};
 const isBodySurface=p=>p.system==='integumentary'&&!(p.id.startsWith('VH_F_')&&p.id!=='VH_F_skin');
 const factory=new Function('T','camera','controls','renderer','ready','latest','atlas','pickers','data','partIsVisible','isBodySurface','bounds','worldBox','hitPoint','hover','dirty',compiled);
 const wheel=factory(T,camera,controls,renderer,true,latest,atlas,pickers,data,partIsVisible,isBodySurface,bounds,new T.Box3(),new T.Vector3(),hover,false);
 const event=(deltaY=-100,extra={})=>({deltaY,deltaMode:0,ctrlKey:false,clientX:rect.left+650,clientY:rect.top+330,prevented:false,stopped:false,preventDefault(){this.prevented=true},stopImmediatePropagation(){this.stopped=true},...extra});
 return {camera,controls,wheel,event,pickers,latest};
};
const screenNear=(a,b,msg)=>assert.ok(Math.hypot(a.x-b.x,a.y-b.y)<1e-9,msg);
const near=(a,b,msg)=>assert.ok(a.distanceTo(b)<1e-9,msg+`: ${a.distanceTo(b)}`);
let checks=0;
for(const viewOffset of [false,true]){
 const ctx=make([{z:0}],{viewOffset});const point=new T.Vector2(.3,.175),ray=new T.Raycaster();ray.setFromCamera(point,ctx.camera);const hit=ray.intersectObject(ctx.pickers[0])[0];assert.ok(hit);const before=hit.point.clone().project(ctx.camera);const e=ctx.event();ctx.wheel(e);screenNear(hit.point.clone().project(ctx.camera),before,'picked surface stays under pointer');assert(e.prevented&&e.stopped);checks++;
}
{
 const ctx=make();const ray=new T.Raycaster();ray.setFromCamera(new T.Vector2(.3,.175),ctx.camera);const anchor=ray.ray.intersectPlane(new T.Plane(new T.Vector3(0,0,1),0),new T.Vector3());const before=anchor.clone().project(ctx.camera);ctx.wheel(ctx.event(-120));screenNear(anchor.clone().project(ctx.camera),before,'empty-space target plane anchor');checks++;
}
for(const hidden of [{state:{visible:['muscular']}},{state:{breastView:'muscle'}},{rendered:false}]){
 const foreground={id:'VH_F_fat_L',system:'mammary',z:1,...(hidden.rendered===false?{rendered:false}:{})};
 const mixed=make([{z:0},foreground],{state:hidden.state});const baseline=make([{z:0}],{state:hidden.state});mixed.wheel(mixed.event());baseline.wheel(baseline.event());near(mixed.camera.position,baseline.camera.position,'hidden breast cannot anchor zoom');checks++;
}
{
 const ctx=make([{z:0},{z:1,system:'integumentary',id:'FJ-skin'}]),solid=make([{z:0}]);ctx.wheel(ctx.event());solid.wheel(solid.event());near(ctx.camera.position,solid.camera.position,'transparent body skin skipped over solids');const skin=make([{z:1,system:'integumentary',id:'FJ-skin'}]);skin.wheel(skin.event());assert(skin.camera.position.distanceTo(solid.camera.position)>1e-3);checks++;
}
{
 const pixels=make(),lines=make(),pages=make(),pinch=make();pixels.wheel(pixels.event(-80));lines.wheel(lines.event(-5,{deltaMode:1}));pages.wheel(pages.event(-.1,{deltaMode:2}));pinch.wheel(pinch.event(-32,{ctrlKey:true}));for(const ctx of [lines,pages,pinch])near(ctx.camera.position,pixels.camera.position,'wheel delta units/pinch normalization');checks++;
}
{
 const ctx=make([],{controls:{minDistance:3,maxDistance:5}});for(let i=0;i<20;i++)ctx.wheel(ctx.event(-10000));assert(Math.abs(ctx.camera.position.distanceTo(ctx.controls.target)-3)<1e-9);for(let i=0;i<20;i++)ctx.wheel(ctx.event(10000));assert(Math.abs(ctx.camera.position.distanceTo(ctx.controls.target)-5)<1e-9);checks++;
}
{
 const ctx=make([],{controls:{enableZoom:false}});const before=ctx.camera.position.clone(),e=ctx.event();ctx.wheel(e);near(ctx.camera.position,before,'disabled zoom unchanged');assert(!e.prevented&&!e.stopped);ctx.controls.enableZoom=true;for(const d of [0,NaN,Infinity])ctx.wheel(ctx.event(d));near(ctx.camera.position,before,'zero/nonfinite delta ignored');checks++;
}
assert(source.includes("removeEventListener('wheel',wheel,true)"));
console.log(`PASS ${checks} actual wheel-handler scenario groups; cleanup capture flag verified.`);
