#version 430

layout (location=0) uniform vec4 fpar;
layout (location=0) out vec4 color;
in vec2 pp;

#define time fpar.x/1000.
#define e vec3(0.0,.001,0.0)
#define normal(DF,p) normalize(vec3(DF((p)+e.yxx).x,DF((p)+e.xyx).x,DF((p)+e.xxy).x)-DF(p).x)

#define la 831u //
//#define lb 9199u
#define lc 207u //
#define ld 524519u
#define le 463u //
#define lf 271u //
#define lg 751u //
#define lh 819u //
#define li 3276u //
//#define lj 65656u
//#define lk 41219u
#define ll 195u //
//#define lm 24627u
#define ln 49203u //
#define lo 255u //
#define lp 799u //
//#define lq 33023u
#define lr 33567u //
#define ls 1006u //
#define lt 3084u //
#define lu 243u //
#define lv 196626u //
#define lw 36915u //
//#define lx 61440u
#define ly 25600u //
//#define lz 12492u
#define l0 3207u //
#define l1 265216u //
#define l2 2437u //
//#define l3 3460u
//#define l4 3330u
//#define l5 1414u
//#define l6 1415u
#define l7 3076u //
#define l8 3463u //
#define l9 3462u //

// GLOBALES --------------------------------------------------------------------------------

vec2 poscol=vec2(0.);
vec3 poscol2 = vec3(0.);
vec3 ppp, pppp;
float quilombo=68.5;

// FUNCIONES GENERALES --------------------------------------

float rand(vec2 p) {return fract(sin(dot(p,vec2(2132.342,4323.343)))*1325.2158);}

vec2 hash(vec2 p) {
	p = vec2(dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)));
	return 2.0*fract(sin(p)*43758.5453123)-1.;
}

float sphere(vec3 p, vec3 rd, float r){
	float b = dot( -p, rd ), i = b*b - dot(p,p) + r*r;
	return i < 0. ?  -1. : b - sqrt(i);
}

vec2 sphere2(vec3 p, vec3 rd, float r){
	float b = dot( -p, rd ),
	inner = b*b - dot(p,p) + r*r;
	return inner < 0. ?  vec2(-1.) : vec2(b - sqrt(inner), b + sqrt(inner));
}


mat2 rot2D(float a) {
    a=radians(a);
    return mat2(cos(a),sin(a),-sin(a),cos(a));
}

mat3 lookat(vec3 fw,vec3 up){
    fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));
    return mat3(rt,cross(rt,fw),fw);
}

float cubo( vec3 p, vec3 b )
{
  return length(max(vec3(0.),abs(p)-b));
}

float toro( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float rect(vec2 p, vec2 b, float inc) {
    p.x+=p.y*inc;
    vec2 d = abs(p)-b;
    return length(max(d,vec2(0))) + min(max(d.x,d.y),0.0);
}

float tri(vec2 p, vec2 q, float ang) {
    p*=rot2D(ang);
    p.x = abs(p.x);
    vec2 a = p - q*clamp( dot(p,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = p - q*vec2( clamp( p.x/q.x, 0.0, 1.0 ), 1.0 );
    float s = -sign( q.y );
    vec2 d = min( vec2( dot(a,a), s*(p.x*q.y-p.y*q.x) ),
                  vec2( dot(b,b), s*(p.y-q.y)  ));
    return -sqrt(d.x)*sign(d.y);
}



// LOGO PVM -----------------------------------------

float logo(vec2 uv) {
    uv.x*=1.15;
    float d=rect(uv,vec2(.045,.25),-.5);
	uv+=vec2(.25,.01);
    d=min(d,rect(uv,vec2(.045,.24),.5));
	uv+=vec2(.265,-.04);
    d=min(d,rect(uv,vec2(.045,.2),-.55));
	uv-=vec2(.73,.06);
    d=min(d,rect(uv,vec2(.045,.16),.4));
	uv+=vec2(-.105,.074);
    d=min(d,rect(uv,vec2(.045,.085),-.45));
	uv+=vec2(-.105,.045);
    d=min(d,rect(uv,vec2(.045,.13),.45));
	uv+=vec2(-.25,.1);
    d=min(d,rect(uv,vec2(.18,.03),.0));
	uv+=vec2(1.32,-.18);
    d=min(d,rect(uv+vec2(0.0,.03),vec2(.35,.03),.0));
	uv+=vec2(-.5165,.4);
    d=min(d,tri(uv,vec2(.09,.185),0.));
	uv-=vec2(.492,.56);
    d=min(d,tri(uv+vec2(0.,.0),vec2(.05,.14),180.)-.01);
    uv+=vec2(.225,-.17);
    d=min(d,tri(uv,vec2(.063,.145),180.));
	uv+=vec2(-.142,.555);
    d=min(d,tri(uv,vec2(.063,.145),0.));
	uv+=vec2(.985,.075);
    vec2 uvd=uv-vec2(uv.y,0.);
    d=min(d,tri(uvd-vec2(0.003,0.022),vec2(.04,.05),0.));
	uv-=vec2(.16,.63);
    uvd=uv+vec2(uv.y*.4,0.);
    d=min(d,tri(uvd+vec2(.03,0.),vec2(.07,.23),-145.));
   	uvd=uv+vec2(.465,.33);
    uvd*=rot2D(27.);
    uvd-=vec2(uvd.y*.5*sign(uvd.x),0.);
    d=min(d,rect(uvd,vec2(.08,.03),.0));
   	uvd=uv+vec2(-1.43,.534);
    uvd*=rot2D(206.);
    uvd-=vec2(uvd.y*.5*sign(uvd.x),0.);
    d=min(d,rect(uvd,vec2(.08,.03),.0));
    float s=pow(abs(d)+.9,10.);
	uvd=uv+vec2(-.28,.36);
    uvd*=rot2D(50.);
    d=max(d,-rect(uvd,vec2(.1,.025),.0));
    return d;
}



// FUNCIONES DE DISTANCIA --------------------------------------------------------------------------------


vec2 de1(vec3 p) {
    p.z += sin(p.y*.6+time*1.5)*3.;
    p.x += cos(p.y*.12)*(.5+sin(p.y*.2)*5.);
    float i=0.;
    float r = length(p);
    ppp = p;
    ppp.y-=time;
    poscol=p.xy;
    pppp = p;
    pppp.xz=rot2D(pppp.y*100.)*pppp.xz;
    pppp=abs(pppp);
    pppp.y+=time*4.;
    pppp.y = abs(.5-fract(pppp.y*.2));
    ppp.y = abs(.5-fract(ppp.y));
	float esf = toro(ppp,vec2(1.8,0.4));
    float cub = cubo(pppp-vec3(2.,0.,0.), vec3(.3));
    cub = min(cub,cubo(pppp-vec3(0.,0.,2.), vec3(.3)));
    float d = esf;
    d = min(d,cub*.3);
    if (abs(esf-d)<.001) i=1.;
    if (abs(cub-d)<.021) i=2.;
    d = max(d,abs(p.y)-100.);
    return vec2(d*.4,i);
}

vec2 de2(vec3 p) {
    p.z += sin(p.y*.3)*2.;
    p.x += cos(p.y*.3)*2.;
    float i=0.;
    float r = length(p);
    ppp = p;
    ppp.y-=time*2.;
    poscol2=p;
    pppp = p;
    pppp.xz=rot2D(pppp.y*20.)*pppp.xz;
    pppp=abs(pppp);
    pppp.y-=time*3.;
    pppp.y = abs(.5-fract(pppp.y*.5));
    ppp.y = abs(.5-fract(ppp.y));
    float esf = toro(ppp,vec2(2.+sin(p.y)*.5,0.1));
    float cub = cubo(pppp-vec3(3.,0.,0.), vec3(.4));
    cub = min(cub,cubo(pppp-vec3(0.,0.,3.), vec3(.4)));
    float d = esf;
    d = min(d,cub*.6);
    if (abs(esf-d)<.001) i=1.;
    if (abs(cub-d)<.001) i=2.;
    return vec2(d*.6,i);
}

vec2 de3(vec3 p) {
    p.z += sin(p.y*1.8)*2.;
    p.x += cos(p.y*1.12)*(.5+sin(p.y*.2)*8.); // # composicion sarasa 3
    float i=0.;
    vec3 prot = p;
    float r = length(p);
    vec3 p3 = p;
    p3.y-=time;
    poscol=p.xy;
    vec3 p4 = p;
    p4.xz=rot2D(p4.y*10.)*p4.xz; // # composisicion sarasa 2
    p4=abs(p4);
    p4.y = abs(.5-fract(p4.y)); // # composicion sarasa
    p3.y = abs(.5-fract(p3.y));
	float esf = toro(p3,vec2(25.8,0.5))*.1; // ver aca # tamanio
    vec3 p2 = p;
    float cub = cubo(p4-vec3(2.,5.,0.), vec3(.3));
    cub = min(cub,cubo(p4-vec3(0.,0.,2.), vec3(.3)))-.5-length(sin(p*8.+p.z*15.))*.05;
    float d = esf;
    d = min(d,cub*.1);
    if (abs(esf-d)<.1) i=1.;
    if (abs(cub-d)<.1) i=2.;
    return vec2(d*.8,i);
}


// ILUMINACION GENERAL --------------------------------------------------------------------------------


vec3 light(vec3 p, vec3 dir, vec3 n, vec3 col, vec3 luzdir, float spec) {	
    return col*(max(0.,dot(dir,-n))*.5+max(0.,dot(luzdir,-n))+max(0.,dot(-luzdir,-n))+0.1)+pow(max(0.,dot(reflect(dir,-n),luzdir))+max(0.,dot(reflect(dir,-n),-luzdir)),5.)*spec;
}

// SHADINGS --------------------------------------------------------------------------------

vec3 shade1(vec3 p, vec3 dir, float i) {
	vec3 n = normal(de1,p);
    vec3 col = vec3(1.-step(.001,abs(i-1.)))*.3;
    col = col*vec3(cos(poscol*5.),0.);
    return light(p, dir, n, col, vec3(0.,0.,1.),1.);
}

vec3 shade2(vec3 p, vec3 dir, float i) {
	vec3 n = normal(de2,p);
    float c = 1.-step(.001,abs(i-1.));
    vec3 colrect = vec3(0.3,cos(p.x)*.3,sin(p.y*.5)*.2);
    vec3 colaros = vec3(abs(sin(poscol2*.5)))*.3;
    vec3 col = mix(colrect,colaros,c);
    return light(p, dir, n, col, vec3(0.,1.,0.),1.);
}

vec3 shade3(vec3 p, vec3 dir, float i) {
	vec3 n = normal(de3,p);
    vec3 col = vec3(1.-step(.001,abs(i-1.)))*.3;
    col *= vec3(sin(poscol*2.),1.);
    col += vec3(1,.7,.4)*.3;
    return light(p, dir, n, col, vec3(0.,0.,1.),.5);
}

// RAYMARCHING CON BACKGROUNDS --------------------------------------------------------------------------------


vec3 march1(vec3 from, vec3 dir) {
    vec3 col=vec3(0.);
    float totdist,glow=0.;
    vec2 d;
    vec3 p;
    for (int i=0; i<150; i++) {
        p=from+totdist*dir;
        d=de1(p);
        if (d.x<.005 || totdist>50.) break;
        totdist+=d.x;
        glow+=1.;
    }
    float l = dot(dir.xy,vec2(0.,1.));
    float luzglow = pow(abs(l),10.)*4.;
    vec3 backcol=vec3(10.,.2,.2)*(mod(pow(abs(l),1.5),.06)+luzglow);
    if (d.x<.01) {
        col=shade1(p-.01*dir,dir,d.y);
        col = mix(col,luzglow*vec3(1.,.8,.7),min(1.,totdist/500.));
    } else {
        col=backcol;
    }
    col +=pow(abs(glow/250.),1.5)*(.7+vec3(1.3,1.3,0.6));
    return col;
}

vec3 march2(vec3 from, vec3 dir) {
    vec3 col=vec3(0.);
    float totdist=0.;
    float st=0.;
    float glow=0.;
    vec2 d;
    vec3 p;
    for (int i=0; i<80; i++) {
        p=from+totdist*dir;
        d=de2(p);
        if (d.x<.005 || totdist>80.) break;
        totdist+=d.x;
        glow+=1.;
    }
    float l = abs(dot(dir.xy,vec2(0.,1.)));
    float luzglow = pow(abs(l),10.)*1.5;
    vec3 colback = vec3(1.5,1.,0.5)*(1.2+cos(dir*25.));
    vec3 colglow = vec3(2.,.7,.3);
    vec3 backcol=colback*(mod(pow(abs(l),3.),.05)+luzglow);
    if (d.x<.1) {
        col=shade2(p-.005*dir,dir,d.y);
        col = mix(col,luzglow*colglow,min(1.,totdist/500.));
    } else {
        col=backcol;
    }
    col +=pow(abs(glow/80.),3.)*colglow;
    return col;
}

vec3 march3(vec3 from, vec3 dir) {
    vec3 col=vec3(0.);
    float totdist=0.;
    float st=0.;
    float glow=0.;
    vec2 d;
    vec3 p;
    for (int i=0; i<220; i++) {
        p=from+totdist*dir;
        d=de3(p);
        if (d.x<0.005 || totdist>100.) break;
        totdist+=d.x;
		glow+=1.;
    }
    float l = dot(dir.xy,vec2(0.,1.));
    float luzglow = pow(abs(l),10.)*1.5;
    vec3 backcol=vec3(luzglow);
	float id = d.y;
    if (d.x<0.1) {
        col=shade3(p-0.005*dir,dir,d.y);
        col = mix(col,luzglow*vec3(1.,.9,.8),min(1.,totdist/70.));
    } else {
        col=backcol;
    }
	//col +=pow(glow/350.,3.)*(1.-step(1.5,id))*vec3(1.,.9,.8);
	col +=pow(abs(glow/400.),2.)*(1.+vec3(1.3,0,0.))*(1.-step(1.5,id));
    return col;
}


// PATHS DE CAMARAS --------------------------------------------------------------------------------


vec3 path(float t){
    t*=.25;
    vec3 p=vec3(cos(t*2.)*4.,-cos(t)*15.,sin(t)*15.);
    return p;
}

vec3 path2(float t){
    t*=.25;
    vec3 p=vec3(cos(t*2.)*3.,-cos(t)*25.,4.-sin(t)*10.);
    return p;
}

vec3 path3(float t){
    t*=.1;
    vec3 p=vec3(5.,-cos(t)*100.,5.); // este es el zoom #
    return p;
}

// SECUENCIAS RAYMARCHING --------------------------------------------------------------------------------

vec3 raymarchers(in vec2 uv, float seq) {
	vec3 col = vec3(0.);
	vec3 dir;
	vec3 from;
	vec3 dircam;
	float ti=time*(1.+step(53.5,time)+step(quilombo,time));
    if (seq<1.) {
		dir = normalize(vec3(uv,.7));
		from = path(ti);
		dircam = normalize(-from);
		dir = lookat(dircam, vec3(1.,1.,0.)) * dir;
		col = march1(from,dir);
	}
    if (abs(seq-2.)<1.) {
		dir = normalize(vec3(uv,1.3));
		from = path2(ti);
		dircam = normalize(-from);
		dir = lookat(dircam, vec3(1.,1.,0.)) * dir;
		col = march2(from,dir);
	}
    if (abs(seq-3.)<1.) {
		uv*=rot2D(time*50.);
		dir = normalize(vec3(uv,.8));
		from = path3(ti);
		dircam = normalize(-from);
		dir = lookat(dircam, vec3(1.,1.,0.)) * dir;
		col = march3(from,dir);
	}
	return col;
}

// KALISET VOLUMETRICO --------------------------------------------------------------------------------------------

vec3 formula(vec3 p) {
	vec3 pos = p;
    p*=1.8;
    float m=100.;
    for (int i=0; i<5; i++) {
    	p=abs(p*3.)/dot(p,p)-1.;
        m=min(max(abs(p.x),abs(p.y)),m);
    }
    m=pow(max(0.,1.-m),5.);
	vec3 col=normalize(pos);
    float ti=time*20.;
    col.xy*=mat2(cos(ti),-sin(ti),sin(ti),cos(ti));
    return abs(col)*m*m*2.+dot(p,p)*.0025;
}

vec3 kalisetVolumetrico(vec2 uv)
{
    float c,s, t=time*.15-2.;
    float distorig=2.+sin(t*3.+56.)*1.3;    
    vec3 ro = -vec3(vec2(sin(t*2.))*.15, distorig),
	     rd =normalize(vec3(uv,3.)),
	     v = vec3(0), p;
    vec2 tt;
 
    mat2 rot;
    c=cos(-t),s=sin(-t);
    rot = mat2(c,-s,s,c);    
    ro.xz*=rot;
    rd.xz*=rot;
    c=cos(-t*1.5),s=sin(-t*1.5);
    rot = mat2(c,-s,s,c);    
    ro.yz*=rot;
    rd.yz*=rot;

    c=cos(t),s=sin(t);
    rot = mat2(c,-s,s,c);    
    ro.xy*=rot;
    rd.xy*=rot;

    float st=.015*100./40.;
	for (float i=1.; i<40.; i++) {
		float d=i*st;
        tt = sphere2(ro, rd, d);
		p = ro+rd*tt.x;
		v+=formula(p)*step(0.,tt.x)*smoothstep(0.,distorig,distorig-tt.x+.4);
		p = ro+rd*tt.y;
		v+=formula(p)*step(0.,tt.x)*smoothstep(0.,distorig,distorig-tt.y+.4);
	}

	return v*.4;
}



// BACKGROUND FLASHERO --------------------------------------------------------------------------

vec3 kset2(vec3 p) {
    vec3 pos=p;
    //p*=1.5;
    vec3 col=vec3(0.);
    float ml=1000.,mc=1000., l, mi;
	for (float i=0.; i<5.; i++) {
		p = abs(p)/dot(p,p) - .5;
        l=min(abs(p.x),abs(p.y));
        if (ml>l) ml=l, mc=length(p),col=normalize(abs(p));
	}
    float c=pow(max(0.,1.0-ml), 80.);
    col=c*abs(normalize(3.+col));
    col*=1.5-pow(mod(l-pos.y+time*.15-.5,1.),.5);
    col+=pow(max(0.,1.-mc),7.);
    return col*min(1.,time);
}

vec3 backgroundFlashero(vec2 uv, float ang) {
	uv*=.6;
    float t=time*.05;
    vec2 mo = vec2(.0,-.25);
    vec3 ro = -vec3(mo, 2.5+cos(time*.2+2.)*.5),
	     rd = normalize(vec3(uv,3.)),
	     v = vec3(0.);
   
	mat2 rot1 = rot2D(ang);
	mat2 rot2 = rot2D(time*2.23+3.);

	for (float i=20.; i<70.; i++) {
		float tt = sphere(ro, rd, i*.02);
		vec3 p = ro+rd*tt;
		p.xz *= rot1;
		p.yz *= rot2;
		vec3 k=kset2(p)*step(0.,tt)*(1.-exp(-.03*(i-20.)));
        //v = .9 * v + .3*k;
        v=max(v,k);

    }

	return v*v;
}

// INTRO -----------------------------------------------------------------------------

vec3 intro(vec2 uv, vec3 background) {
    vec2 p = uv;
	p-=vec2(0.1,0.);
    float d=logo(p)-.02;
	for (int i=0; i< 20; i++) {
        vec2 ppp = hash(vec2(float(i), float(i*2)));
        float x = length(p*1000. - ppp);
        if (x < 10.1) d*=x;
	}

	vec3 col = background * (1.0 - exp(-4.0*(pow(abs(d), 1.2))));
	// Affeccion de lo blanco por el fondo
    col = mix(col, vec3(1.,.5,.3), 1.0 - smoothstep(0.0,0.025,abs(d)));
    // Engorde del nodo
    col += vec3(2.,.5,.3)*(mix(col, vec3(0.5), .9 - smoothstep(0.0,0.1,abs(d-background.r/3.5))));
    // 3D??
    //col = max(col,.7*vec3(1.,1.,1.2)*(1.-smoothstep(-.1,.0,logo(p+vec2(.04,-.04)))));
    return col*min(1.,time);;
}

// LETRAS ------------------------------------------------------------------------------------------------

uint frase1[]=uint[](lh,li,0u,ln,li,lc,le,0u,lp,le,lo,lp,ll,le,0u,la,lt,0u,lr,le,lv,li,ls,li,lo,ln,999u);
uint frase2[]=uint[](lt,lh,li,ls,0u,lw,la,ls,0u,ls,lp,la,lc,le,0u,lp,lo,lo,0u,0u,999u);
uint frase3[]=uint[](lc,0u,lu,0u,la,lg,la,li,ln,0u,la,lt,0u,lf,ll,la,ls,lh,lp,la,lr,lt,ly,0u,l2,l0,l1,l9,999u);
uint frase4[]=uint[](ls,le,lp,lt,0u,l2,l7,0u,l2,l8,0u,l2,l9,0u,0u,la,lr,lg,le,ln,lt,li,ln,la,999u);

float linea(vec2 uv, vec4 vert){
    uv-=vert.xy;
    vec2 pos=vert.zw;
    return length(uv-pos*clamp(dot(uv,pos)/dot(pos,pos),0.,1.));
}

float letra(vec2 uv, uint letra, float width) {
    uv*=1.+length(uv-vec2(.5,.7))*.7;
    vec4 d = vec4(0.5,-0.5,1.,0.);
	float s = 10.;
    s = (letra & 1u     ) > 0u ? min(s,linea(uv, d.wwwx)) : s; // vertical -1, -1
	s = (letra & 2u     ) > 0u ? min(s,linea(uv, d.wxwx)) : s; // vertical -1,  1
	s = (letra & 4u     ) > 0u ? min(s,linea(uv, d.wzxw)) : s; // horizont -1,  1
	s = (letra & 8u     ) > 0u ? min(s,linea(uv, d.xzxw)) : s; // horizont  1,  1
 	s = (letra & 16u    ) > 0u ? min(s,linea(uv, d.zxwx)) : s; // vertical  1,  1
 	s = (letra & 32u    ) > 0u ? min(s,linea(uv, d.zwwx)) : s; // vertical  1, -1
 	s = (letra & 64u    ) > 0u ? min(s,linea(uv, d.xwxw)) : s; // horizont  1, -1
 	s = (letra & 128u   ) > 0u ? min(s,linea(uv, d.wwxw)) : s; // horizont -1, -1
 	s = (letra & 256u   ) > 0u ? min(s,linea(uv, d.wxxw)) : s; // horizont -1,  0
 	s = (letra & 512u   ) > 0u ? min(s,linea(uv, d.xxxw)) : s; // horizont  1,  0
 	s = (letra & 1024u  ) > 0u ? min(s,linea(uv, d.xwwx)) : s; // vertical  0, -1
 	s = (letra & 2048u  ) > 0u ? min(s,linea(uv, d.xxwx)) : s; // vertical  0,  1
 	s = (letra & 4096u  ) > 0u ? min(s,linea(uv, d.wwxx)) : s; // diagonal -1, -1
 	s = (letra & 8192u  ) > 0u ? min(s,linea(uv, d.xxxx)) : s; // diagonal  1,  1
 	s = (letra & 16384u ) > 0u ? min(s,linea(uv, d.wzxy)) : s; // diagonal -1,  1
 	s = (letra & 32768u ) > 0u ? min(s,linea(uv, d.xxxy)) : s; // diagonal  1, -1
 	s = (letra & 65536u ) > 0u ? min(s,linea(uv, d.wxxy)) : s; // diagona2 -1, -1
 	s = (letra & 131072u) > 0u ? min(s,linea(uv, d.xwxx)) : s; // diagona2  1, -1
 	s = (letra & 262144u) > 0u ? min(s,linea(uv, d.wxxx)) : s; // diagona2 -1,  1
 	s = (letra & 524288u) > 0u ? min(s,linea(uv, d.xzxy)) : s; // diagona2  1,  1
    return s-width;
}

// OUTRO ---------------------------------------------------------------------

vec3 textos(vec2 uv, float id, float scroll, vec3 background) {
	vec2 p = uv;
    if (id<1.) p.x+=1.4;
	if (id>1. && id<2.) p.x+=1.1;
	if (id>2. && id<3.) p.x+=1.5;
	if (id>3.) p.x+=1.3;
	p.x-=scroll*5.;
    p*=10.;
    vec3 c=vec3(0.);
	uint let;
    for (int i=0; i<30; i++) {
        p.y+=sin(float(i)*2.+time*2.)*.2;
    	if (id<4.) let = frase4[i];
    	if (id<3.) let = frase3[i];
    	if (id<2.) let = frase2[i];
    	if (id<1.) let = frase1[i];
        if (let == 999u) break;
        c+=pow(max(0., 1.-letra(p, let, .1)),13.);
        p.x-=1.1;
    }
    return background+c*min(1.,time*.3);
}

// --------------------------------------------------------------------------------


void main() {
	vec2 uv=pp;
	vec3 c=vec3(0.);
	uv.x*=1.77;
	float seq=mod(time*0.258333333+.25,3.);
	float id = mod((time-84.)*.18,4.);
	float scroll=smoothstep(.8,1.,fract(id));
	float ti=time*(1.+step(quilombo,time));
	float fr=abs(.5-fract((ti-53.)*.2583333333*4.));
	float flash=step(53.,time)*pow(smoothstep(0.3,0.5,fr),2.-step(quilombo,time)*1.5);
    float ang=90.*floor(id)+scroll*90.;
	vec3 background = backgroundFlashero(uv, ang*step(10.,time));
	if (time<7.5) 
		c=intro(uv, background);		
	if ((time>7.7 && time<45.) || (time>53.5 && time<83.5)) {
		uv=pow(abs(uv),1.+vec2(flash))*sign(uv);
		c=raymarchers(uv, seq)*(1.+flash*3.)+smoothstep(78.,83.5,time)*mod(time,.1)/.1;
	}
	if (time > 45.5 && time < 53.5) 
		c=kalisetVolumetrico(uv);
	if (time >84.) 
		c=textos(uv,id,scroll,background);
	color=vec4(c,1.);
}
