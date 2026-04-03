var M=Object.defineProperty;var S=(l,r,e)=>r in l?M(l,r,{enumerable:!0,configurable:!0,writable:!0,value:e}):l[r]=e;var g=(l,r,e)=>S(l,typeof r!="symbol"?r+"":r,e);const U=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`;const N=`
precision highp float;
precision highp int;

uniform sampler2D u_source;
uniform vec2 u_size;
uniform vec2 u_direction;
uniform int u_sourceMode;
uniform float u_kernel[33];

float unpackByte(float value) {
  return floor(value * 255.0 + 0.5);
}

float unpackRawLA(vec4 sample) {
  float lo = unpackByte(sample.r);
  float hi = unpackByte(sample.a);
  return (lo + hi * 256.0) / 65535.0;
}

float unpackRG16(vec4 sample) {
  float lo = unpackByte(sample.r);
  float hi = unpackByte(sample.g);
  return (lo + hi * 256.0) / 65535.0;
}

vec2 pixelUv(float x, float y) {
  vec2 clamped = clamp(vec2(x, y), vec2(0.0), u_size - 1.0);
  return (clamped + 0.5) / u_size;
}

float sampleSource(float x, float y) {
  vec4 sample = texture2D(u_source, pixelUv(x, y));
  return u_sourceMode == 0 ? unpackRawLA(sample) : unpackRG16(sample);
}

vec2 packRG16(float value) {
  float scaled = floor(clamp(value, 0.0, 1.0) * 65535.0 + 0.5);
  float lo = mod(scaled, 256.0);
  float hi = floor(scaled / 256.0);
  return vec2(lo, hi) / 255.0;
}

void main() {
  float x = floor(gl_FragCoord.x - 0.5);
  float y = floor(gl_FragCoord.y - 0.5);
  float sum = 0.0;
  for (int i = 0; i < 33; i++) {
    float offset = float(i - 16);
    sum += sampleSource(x + u_direction.x * offset, y + u_direction.y * offset) * u_kernel[i];
  }
  vec2 packed = packRG16(sum);
  gl_FragColor = vec4(packed.x, packed.y, 0.0, 1.0);
}
`,P=`
precision highp float;
precision highp int;

uniform sampler2D u_blurred;
uniform vec2 u_size;
uniform float u_maxSearchRadius;

const float TWO_PI = 6.28318530718;
const float MIN_GAP = 0.523598775598;
const float ANGLE_STEP = TWO_PI / float(36);

float unpackByte(float value) {
  return floor(value * 255.0 + 0.5);
}

float unpackRG16(vec4 sample) {
  float lo = unpackByte(sample.r);
  float hi = unpackByte(sample.g);
  return (lo + hi * 256.0) / 65535.0;
}

vec2 pixelUv(float x, float y) {
  vec2 clamped = clamp(vec2(x, y), vec2(0.0), u_size - 1.0);
  return (clamped + 0.5) / u_size;
}

float sampleBlurred(float x, float y) {
  return unpackRG16(texture2D(u_blurred, pixelUv(x, y)));
}

float circularDistance(float a, float b) {
  float diff = abs(a - b);
  return min(diff, TWO_PI - diff);
}

vec2 packDistance(float value) {
  float scaled = floor(clamp(value, 0.0, 65535.0) + 0.5);
  float lo = mod(scaled, 256.0);
  float hi = floor(scaled / 256.0);
  return vec2(lo, hi) / 255.0;
}

void main() {
  float x = floor(gl_FragCoord.x - 0.5);
  float y = floor(gl_FragCoord.y - 0.5);
  float center = sampleBlurred(x, y);
  float safeEdge = floor(min(min(x, y), min(u_size.x - 1.0 - x, u_size.y - 1.0 - y)));
  safeEdge = min(safeEdge, u_maxSearchRadius);

  int bestRadius = 0;
  for (int radius = 1; radius <= 256; radius++) {
    if (float(radius) > safeEdge + 0.5) break;

    int selectedCount = 0;
    float a0 = 0.0;
    float a1 = 0.0;
    float a2 = 0.0;
    float a3 = 0.0;
    float a4 = 0.0;
    float a5 = 0.0;

    for (int angleIndex = 0; angleIndex < 36; angleIndex++) {
      float angle = float(angleIndex) * ANGLE_STEP;
      float sampleValue = sampleBlurred(
        x + cos(angle) * float(radius),
        y + sin(angle) * float(radius)
      );
      if (sampleValue + (1.0 / 65535.0) < center) continue;

      bool ok = true;
      if (selectedCount > 0 && circularDistance(angle, a0) < MIN_GAP) ok = false;
      if (selectedCount > 1 && circularDistance(angle, a1) < MIN_GAP) ok = false;
      if (selectedCount > 2 && circularDistance(angle, a2) < MIN_GAP) ok = false;
      if (selectedCount > 3 && circularDistance(angle, a3) < MIN_GAP) ok = false;
      if (selectedCount > 4 && circularDistance(angle, a4) < MIN_GAP) ok = false;
      if (selectedCount > 5 && circularDistance(angle, a5) < MIN_GAP) ok = false;
      if (!ok) continue;

      if (selectedCount == 0) a0 = angle;
      else if (selectedCount == 1) a1 = angle;
      else if (selectedCount == 2) a2 = angle;
      else if (selectedCount == 3) a3 = angle;
      else if (selectedCount == 4) a4 = angle;
      else if (selectedCount == 5) a5 = angle;
      selectedCount += 1;

      if (selectedCount >= 6) {
        bestRadius = radius;
        break;
      }
    }
  }

  vec2 packed = packDistance(float(bestRadius));
  gl_FragColor = vec4(packed.x, packed.y, 0.0, 1.0);
}
`,I=l=>Math.min(16,Math.max(1,Math.ceil(Math.max(.5,l)*3))),v=l=>{const r=I(l),e=new Float32Array(33),a=Math.max(.5,l);let s=0;for(let n=-r;n<=r;n++){const t=Math.exp(-(n*n)/(2*a*a));e[n+16]=t,s+=t}for(let n=0;n<e.length;n++)e[n]/=Math.max(s,1e-12);return e};class y{constructor(){g(this,"canvas",null);g(this,"gl",null);g(this,"blurProgram",null);g(this,"structureProgram",null);g(this,"positionBuffer",null);g(this,"blurUniforms",null);g(this,"structureUniforms",null);g(this,"resources",null);g(this,"initialized",!1)}async init(){if(this.initialized&&this.gl&&this.blurProgram&&this.structureProgram)return!0;const r=this.createCanvas();if(!r)return!1;const e=r.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return!1;const a=this.compileShader(e,e.VERTEX_SHADER,U),s=this.compileShader(e,e.FRAGMENT_SHADER,N),n=this.compileShader(e,e.FRAGMENT_SHADER,P);if(!a||!s||!n)return a&&e.deleteShader(a),s&&e.deleteShader(s),n&&e.deleteShader(n),!1;const t=this.createProgram(e,a,s),o=this.createProgram(e,a,n);if(e.deleteShader(a),e.deleteShader(s),e.deleteShader(n),!t||!o)return t&&e.deleteProgram(t),o&&e.deleteProgram(o),!1;const u=e.createBuffer();return u?(e.bindBuffer(e.ARRAY_BUFFER,u),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=r,this.gl=e,this.blurProgram=t,this.structureProgram=o,this.positionBuffer=u,this.blurUniforms={source:e.getUniformLocation(t,"u_source"),size:e.getUniformLocation(t,"u_size"),direction:e.getUniformLocation(t,"u_direction"),sourceMode:e.getUniformLocation(t,"u_sourceMode"),kernel:e.getUniformLocation(t,"u_kernel[0]")},this.structureUniforms={blurred:e.getUniformLocation(o,"u_blurred"),size:e.getUniformLocation(o,"u_size"),maxSearchRadius:e.getUniformLocation(o,"u_maxSearchRadius")},this.initialized=!0,!0):(e.deleteProgram(t),e.deleteProgram(o),!1)}async computeDistances(r,e,a,s,n){if(s>256||!this.initialized&&!await this.init())return null;const t=this.gl,o=this.blurProgram,u=this.structureProgram,i=this.blurUniforms,f=this.structureUniforms;if(!t||!o||!u||!i||!f||!this.positionBuffer||!this.canvas)return null;const c=this.ensureResources(e,a);if(!c)return null;const E=this.packRawToLa16(r),m=v(n);this.canvas.width=e,this.canvas.height=a,t.viewport(0,0,e,a),t.disable(t.BLEND),t.pixelStorei(t.UNPACK_ALIGNMENT,1),t.pixelStorei(t.PACK_ALIGNMENT,1),t.bindBuffer(t.ARRAY_BUFFER,this.positionBuffer),t.bindTexture(t.TEXTURE_2D,c.rawTexture),t.texImage2D(t.TEXTURE_2D,0,t.LUMINANCE_ALPHA,e,a,0,t.LUMINANCE_ALPHA,t.UNSIGNED_BYTE,E),t.useProgram(o),t.enableVertexAttribArray(0),t.vertexAttribPointer(0,2,t.FLOAT,!1,0,0),t.uniform2f(i.size,e,a),t.uniform1fv(i.kernel,m),t.bindFramebuffer(t.FRAMEBUFFER,c.pingFramebuffer),t.activeTexture(t.TEXTURE0),t.bindTexture(t.TEXTURE_2D,c.rawTexture),t.uniform1i(i.source,0),t.uniform2f(i.direction,1,0),t.uniform1i(i.sourceMode,0),t.drawArrays(t.TRIANGLES,0,6),t.bindFramebuffer(t.FRAMEBUFFER,c.pongFramebuffer),t.bindTexture(t.TEXTURE_2D,c.pingTexture),t.uniform2f(i.direction,0,1),t.uniform1i(i.sourceMode,1),t.drawArrays(t.TRIANGLES,0,6),t.disableVertexAttribArray(0),t.useProgram(u),t.enableVertexAttribArray(0),t.vertexAttribPointer(0,2,t.FLOAT,!1,0,0),t.bindFramebuffer(t.FRAMEBUFFER,c.outputFramebuffer),t.activeTexture(t.TEXTURE0),t.bindTexture(t.TEXTURE_2D,c.pongTexture),t.uniform1i(f.blurred,0),t.uniform2f(f.size,e,a),t.uniform1f(f.maxSearchRadius,Math.max(1,Math.round(s))),t.drawArrays(t.TRIANGLES,0,6);const p=new Uint8Array(e*a*4);t.readPixels(0,0,e,a,t.RGBA,t.UNSIGNED_BYTE,p),t.disableVertexAttribArray(0),t.bindFramebuffer(t.FRAMEBUFFER,null),t.bindBuffer(t.ARRAY_BUFFER,null),t.bindTexture(t.TEXTURE_2D,null);const d=new Uint16Array(e*a);for(let A=0,R=0;A<d.length;A++,R+=4)d[A]=p[R]+(p[R+1]<<8);return d}packRawToLa16(r){let e=Number.POSITIVE_INFINITY,a=Number.NEGATIVE_INFINITY;for(let t=0;t<r.length;t++){const o=r[t];Number.isFinite(o)&&(e=Math.min(e,o),a=Math.max(a,o))}(!Number.isFinite(e)||!Number.isFinite(a)||a<=e)&&(e=0,a=1);const s=65535/Math.max(a-e,1e-12),n=new Uint8Array(r.length*2);for(let t=0;t<r.length;t++){const o=Number.isFinite(r[t])?r[t]:e,u=Math.max(0,Math.min(65535,Math.round((o-e)*s))),i=t*2;n[i]=u&255,n[i+1]=u>>8&255}return n}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(r,e){const a=this.gl;if(!a)return null;if(this.resources&&this.resources.width===r&&this.resources.height===e)return this.resources;this.disposeResources();const s=this.createTexture(a.LUMINANCE_ALPHA,r,e,a.LUMINANCE_ALPHA,a.UNSIGNED_BYTE,null),n=this.createTexture(a.RGBA,r,e,a.RGBA,a.UNSIGNED_BYTE,null),t=this.createTexture(a.RGBA,r,e,a.RGBA,a.UNSIGNED_BYTE,null),o=this.createTexture(a.RGBA,r,e,a.RGBA,a.UNSIGNED_BYTE,null),u=this.createFramebuffer(n),i=this.createFramebuffer(t),f=this.createFramebuffer(o);return!s||!n||!t||!o||!u||!i||!f?(s&&a.deleteTexture(s),n&&a.deleteTexture(n),t&&a.deleteTexture(t),o&&a.deleteTexture(o),u&&a.deleteFramebuffer(u),i&&a.deleteFramebuffer(i),f&&a.deleteFramebuffer(f),null):(this.resources={width:r,height:e,rawTexture:s,pingTexture:n,pongTexture:t,outputTexture:o,pingFramebuffer:u,pongFramebuffer:i,outputFramebuffer:f},this.resources)}createTexture(r,e,a,s,n,t){const o=this.gl;if(!o)return null;const u=o.createTexture();return u?(o.bindTexture(o.TEXTURE_2D,u),o.texParameteri(o.TEXTURE_2D,o.TEXTURE_MIN_FILTER,o.NEAREST),o.texParameteri(o.TEXTURE_2D,o.TEXTURE_MAG_FILTER,o.NEAREST),o.texParameteri(o.TEXTURE_2D,o.TEXTURE_WRAP_S,o.CLAMP_TO_EDGE),o.texParameteri(o.TEXTURE_2D,o.TEXTURE_WRAP_T,o.CLAMP_TO_EDGE),o.texImage2D(o.TEXTURE_2D,0,r,e,a,0,s,n,t),o.bindTexture(o.TEXTURE_2D,null),u):null}createFramebuffer(r){const e=this.gl;if(!e||!r)return null;const a=e.createFramebuffer();if(!a)return null;e.bindFramebuffer(e.FRAMEBUFFER,a),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,r,0);const s=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),s!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(a),null):a}compileShader(r,e,a){const s=r.createShader(e);return s?(r.shaderSource(s,a),r.compileShader(s),r.getShaderParameter(s,r.COMPILE_STATUS)?s:(console.error("[StructureLength WebGL] shader compile failed",r.getShaderInfoLog(s)),r.deleteShader(s),null)):null}createProgram(r,e,a){const s=r.createProgram();return s?(r.attachShader(s,e),r.attachShader(s,a),r.bindAttribLocation(s,0,"a_position"),r.linkProgram(s),r.getProgramParameter(s,r.LINK_STATUS)?s:(console.error("[StructureLength WebGL] program link failed",r.getProgramInfoLog(s)),r.deleteProgram(s),null)):null}disposeResources(){const r=this.gl,e=this.resources;if(!r||!e){this.resources=null;return}r.deleteTexture(e.rawTexture),r.deleteTexture(e.pingTexture),r.deleteTexture(e.pongTexture),r.deleteTexture(e.outputTexture),r.deleteFramebuffer(e.pingFramebuffer),r.deleteFramebuffer(e.pongFramebuffer),r.deleteFramebuffer(e.outputFramebuffer),this.resources=null}}const D=new y,b=(l,r)=>{if(r<=1)return 0;let e=l;const a=r-1;for(;e<0||e>a;)e=e<0?-e:a-(e-a);return e},L=l=>{const r=Math.max(.5,l),e=Math.max(1,Math.ceil(r*3)),a=new Float32Array(e*2+1);let s=0;for(let n=-e;n<=e;n++){const t=Math.exp(-(n*n)/(2*r*r));a[n+e]=t,s+=t}for(let n=0;n<a.length;n++)a[n]/=Math.max(s,1e-12);return a},k=(l,r,e,a)=>{if(a<=1e-6)return new Float32Array(l);const s=L(a),n=s.length-1>>1,t=new Float32Array(l.length),o=new Float32Array(l.length);for(let u=0;u<e;u++){const i=u*r;for(let f=0;f<r;f++){let c=0;for(let E=-n;E<=n;E++){const m=b(f+E,r);c+=l[i+m]*s[E+n]}t[i+f]=c}}for(let u=0;u<e;u++)for(let i=0;i<r;i++){let f=0;for(let c=-n;c<=n;c++){const E=b(u+c,e);f+=t[E*r+i]*s[c+n]}o[u*r+i]=f}return o},C=(l,r,e,a,s)=>{const n=Math.max(0,Math.min(r-1,a)),t=Math.max(0,Math.min(e-1,s)),o=Math.floor(n),u=Math.floor(t),i=Math.min(r-1,o+1),f=Math.min(e-1,u+1),c=n-o,E=t-u,m=l[u*r+o],p=l[u*r+i],d=l[f*r+o],A=l[f*r+i],R=m+(p-m)*c,h=d+(A-d)*c;return R+(h-R)*E},B=(l,r,e)=>{if(l.length<e)return!1;const a=[...l].sort((n,t)=>n-t),s=a.concat(a.map(n=>n+Math.PI*2));for(let n=0;n<a.length;n++){let t=1,o=s[n];for(let u=n+1;u<n+a.length;u++){const i=s[u];if(!(i-o+1e-6<r)&&(t++,o=i,t>=e)){if(a[n]+Math.PI*2-o+1e-6>=r)return!0;break}}}return!1},F=l=>{const r=l.length;if(r<=0)return 0;if(r===1)return l[0];if(r===2)return(l[0]+l[1])*.5;if(r===3){const n=l[1]-l[0],t=l[2]-l[1];return n<t?(l[0]+l[1])*.5:t<n?(l[1]+l[2])*.5:l[1]}const e=Math.ceil(r*.5);let a=0,s=Number.POSITIVE_INFINITY;for(let n=0;n<=r-e;n++){const t=l[n+e-1]-l[n];t<s-1e-9&&(s=t,a=n)}return F(l.slice(a,a+e))},G=(l,r,e,a,s,n)=>{let t=0,o=0,u=0,i=0;const f=[],c=Math.max(1,Math.round(a)),E=new Float32Array(c+1),m=new Float32Array(c+1);for(let d=0;d<m.length;d++)m[d]=d;for(let d=0;d<l.length;d++){const A=l[d],R=d%r,h=Math.floor(d/r);if(Math.floor(Math.min(R,h,r-1-R,e-1-h,c))<1)continue;i++;const _=Math.max(0,Math.min(c,A));E[_]+=1,f.push(_),_>t&&(t=_,o=R,u=h)}f.sort((d,A)=>d-A);const p=F(f);return{backend:n,sigmaPx:s,maxSearchRadiusPx:c,hsmModePx:p,bestDistancePx:t,bestX:o,bestY:u,evaluatedPixelCount:i,histogramXs:m,histogramCounts:E}},X=(l,r,e,a,s)=>{const n=k(l,r,e,s),t=96,o=new Float32Array(t),u=new Float32Array(t),i=new Float32Array(t);for(let m=0;m<t;m++){const p=m/t*Math.PI*2;o[m]=p,u[m]=Math.cos(p),i[m]=Math.sin(p)}const f=Math.PI/6,c=Math.max(1,Math.round(a)),E=new Uint16Array(r*e);for(let m=0;m<e;m++)for(let p=0;p<r;p++){const d=m*r+p,A=n[d];if(!Number.isFinite(A)){E[d]=0;continue}const R=Math.floor(Math.min(p,m,r-1-p,e-1-m,c));if(R<1){E[d]=0;continue}let h=0;const T=[];for(let _=1;_<=R;_++){T.length=0;for(let x=0;x<t;x++)C(n,r,e,p+u[x]*_,m+i[x]*_)+1e-6>=A&&T.push(o[x]);B(T,f,6)&&(h=_)}E[d]=h}return E};self.onmessage=async l=>{const{roiData:r,width:e,height:a,maxSearchRadiusPx:s,sigmaPx:n}=l.data;try{let t="cpu",o=await D.computeDistances(r,e,a,s,n);o?t="webgl":o=X(r,e,a,s,n);const u=G(o,e,a,s,n,t),i={success:!0,result:u};self.postMessage(i,[u.histogramXs.buffer,u.histogramCounts.buffer])}catch(t){const o={success:!1,error:(t==null?void 0:t.message)||String(t)};self.postMessage(o)}};
