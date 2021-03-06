/* File generated with Shader Minifier 1.1.6
 * http://www.ctrl-alt-test.fr
 */
#ifndef FRAGMENT_SHADER_INL_
# define FRAGMENT_SHADER_INL_

const char *fragment_shader_glsl =
 "#version 430\n"
 "layout(location=0)uniform vec4 fpar;"
 "layout(location=0)out vec4 color;"
 "in vec2 pp;\n"
 "#define time fpar.x/1000.\n"
 "#define e vec3(0.0,.001,0.0)\n"
 "#define normal(DF,p)normalize(vec3(DF((p)+e.yxx).x,DF((p)+e.xyx).x,DF((p)+e.xxy).x)-DF(p).x)\n"
 "#define la 831u\n"
 "#define lc 207u\n"
 "#define ld 524519u\n"
 "#define le 463u\n"
 "#define lf 271u\n"
 "#define lg 751u\n"
 "#define lh 819u\n"
 "#define li 3276u\n"
 "#define ll 195u\n"
 "#define ln 49203u\n"
 "#define lo 255u\n"
 "#define lp 799u\n"
 "#define lr 33567u\n"
 "#define ls 1006u\n"
 "#define lt 3084u\n"
 "#define lu 243u\n"
 "#define lv 196626u\n"
 "#define lw 36915u\n"
 "#define ly 25600u\n"
 "#define l0 3207u\n"
 "#define l1 265216u\n"
 "#define l2 2437u\n"
 "#define l7 3076u\n"
 "#define l8 3463u\n"
 "#define l9 3462u\n"
 "vec2 v=vec2(0.);"
 "vec3 f=vec3(0.),l,y;"
 "float a=68.5;"
 "float n(vec2 v)"
 "{"
   "return fract(sin(dot(v,vec2(2132.34,4323.34)))*1325.22);"
 "}"
 "vec2 s(vec2 v)"
 "{"
   "return v=vec2(dot(v,vec2(127.1,311.7)),dot(v,vec2(269.5,183.3))),2.*fract(sin(v)*43758.5)-1.;"
 "}"
 "float n(vec3 v,vec3 l,float y)"
 "{"
   "float m=dot(-v,l),a=m*m-dot(v,v)+y*y;"
   "return a<0.?-1.:m-sqrt(a);"
 "}"
 "vec2 s(vec3 v,vec3 l,float y)"
 "{"
   "float m=dot(-v,l),a=m*m-dot(v,v)+y*y;"
   "return a<0.?vec2(-1.):vec2(m-sqrt(a),m+sqrt(a));"
 "}"
 "mat2 t(float v)"
 "{"
   "return v=radians(v),mat2(cos(v),sin(v),-sin(v),cos(v));"
 "}"
 "mat3 n(vec3 v,vec3 l)"
 "{"
   "v=normalize(v);"
   "vec3 a=normalize(cross(v,normalize(l)));"
   "return mat3(a,cross(a,v),v);"
 "}"
 "float s(vec3 v,vec3 l)"
 "{"
   "return length(max(vec3(0.),abs(v)-l));"
 "}"
 "float t(vec3 v,vec2 l)"
 "{"
   "vec2 y=vec2(length(v.xz)-l.x,v.y);"
   "return length(y)-l.y;"
 "}"
 "float t(vec2 v,vec2 l,float y)"
 "{"
   "v.x+=v.y*y;"
   "vec2 a=abs(v)-l;"
   "return length(max(a,vec2(0)))+min(max(a.x,a.y),0.);"
 "}"
 "float p(vec2 v,vec2 l,float y)"
 "{"
   "v*=t(y);"
   "v.x=abs(v.x);"
   "vec2 a=v-l*clamp(dot(v,l)/dot(l,l),0.,1.),u=v-l*vec2(clamp(v.x/l.x,0.,1.),1.);"
   "float f=-sign(l.y);"
   "vec2 m=min(vec2(dot(a,a),f*(v.x*l.y-v.y*l.x)),vec2(dot(u,u),f*(v.y-l.y)));"
   "return-sqrt(m.x)*sign(m.y);"
 "}"
 "float p(vec2 v)"
 "{"
   "v.x*=1.15;"
   "float y=t(v,vec2(.045,.25),-.5);"
   "v+=vec2(.25,.01);"
   "y=min(y,t(v,vec2(.045,.24),.5));"
   "v+=vec2(.265,-.04);"
   "y=min(y,t(v,vec2(.045,.2),-.55));"
   "v-=vec2(.73,.06);"
   "y=min(y,t(v,vec2(.045,.16),.4));"
   "v+=vec2(-.105,.074);"
   "y=min(y,t(v,vec2(.045,.085),-.45));"
   "v+=vec2(-.105,.045);"
   "y=min(y,t(v,vec2(.045,.13),.45));"
   "v+=vec2(-.25,.1);"
   "y=min(y,t(v,vec2(.18,.03),0.));"
   "v+=vec2(1.32,-.18);"
   "y=min(y,t(v+vec2(0.,.03),vec2(.35,.03),0.));"
   "v+=vec2(-.5165,.4);"
   "y=min(y,p(v,vec2(.09,.185),0.));"
   "v-=vec2(.492,.56);"
   "y=min(y,p(v+vec2(0.,0.),vec2(.05,.14),180.)-.01);"
   "v+=vec2(.225,-.17);"
   "y=min(y,p(v,vec2(.063,.145),180.));"
   "v+=vec2(-.142,.555);"
   "y=min(y,p(v,vec2(.063,.145),0.));"
   "v+=vec2(.985,.075);"
   "vec2 l=v-vec2(v.y,0.);"
   "y=min(y,p(l-vec2(.003,.022),vec2(.04,.05),0.));"
   "v-=vec2(.16,.63);"
   "l=v+vec2(v.y*.4,0.);"
   "y=min(y,p(l+vec2(.03,0.),vec2(.07,.23),-145.));"
   "l=v+vec2(.465,.33);"
   "l*=t(27.);"
   "l-=vec2(l.y*.5*sign(l.x),0.);"
   "y=min(y,t(l,vec2(.08,.03),0.));"
   "l=v+vec2(-1.43,.534);"
   "l*=t(206.);"
   "l-=vec2(l.y*.5*sign(l.x),0.);"
   "y=min(y,t(l,vec2(.08,.03),0.));"
   "float a=pow(abs(y)+.9,10.);"
   "l=v+vec2(-.28,.36);"
   "l*=t(50.);"
   "y=max(y,-t(l,vec2(.1,.025),0.));"
   "return y;"
 "}"
 "vec2 m(vec3 m)"
 "{"
   "m.z+=sin(m.y*.6+time*1.5)*3.;"
   "m.x+=cos(m.y*.12)*(.5+sin(m.y*.2)*5.);"
   "float a=0.,f=length(m);"
   "l=m;"
   "l.y-=time;"
   "v=m.xy;"
   "y=m;"
   "y.xz=t(y.y*100.)*y.xz;"
   "y=abs(y);"
   "y.y+=time*4.;"
   "y.y=abs(.5-fract(y.y*.2));"
   "l.y=abs(.5-fract(l.y));"
   "float p=t(l,vec2(1.8,.4)),i=s(y-vec3(2.,0.,0.),vec3(.3));"
   "i=min(i,s(y-vec3(0.,0.,2.),vec3(.3)));"
   "float u=p;"
   "u=min(u,i*.3);"
   "if(abs(p-u)<.001)"
     "a=1.;"
   "if(abs(i-u)<.021)"
     "a=2.;"
   "u=max(u,abs(m.y)-100.);"
   "return vec2(u*.4,a);"
 "}"
 "vec2 x(vec3 v)"
 "{"
   "v.z+=sin(v.y*.3)*2.;"
   "v.x+=cos(v.y*.3)*2.;"
   "float a=0.,p=length(v);"
   "l=v;"
   "l.y-=time*2.;"
   "f=v;"
   "y=v;"
   "y.xz=t(y.y*20.)*y.xz;"
   "y=abs(y);"
   "y.y-=time*3.;"
   "y.y=abs(.5-fract(y.y*.5));"
   "l.y=abs(.5-fract(l.y));"
   "float m=t(l,vec2(2.+sin(v.y)*.5,.1)),i=s(y-vec3(3.,0.,0.),vec3(.4));"
   "i=min(i,s(y-vec3(0.,0.,3.),vec3(.4)));"
   "float u=m;"
   "u=min(u,i*.6);"
   "if(abs(m-u)<.001)"
     "a=1.;"
   "if(abs(i-u)<.001)"
     "a=2.;"
   "return vec2(u*.6,a);"
 "}"
 "vec2 u(vec3 y)"
 "{"
   "y.z+=sin(y.y*1.8)*2.;"
   "y.x+=cos(y.y*1.12)*(.5+sin(y.y*.2)*8.);"
   "float l=0.;"
   "vec3 a=y;"
   "float m=length(y);"
   "vec3 f=y;"
   "f.y-=time;"
   "v=y.xy;"
   "vec3 i=y;"
   "i.xz=t(i.y*10.)*i.xz;"
   "i=abs(i);"
   "i.y=abs(.5-fract(i.y));"
   "f.y=abs(.5-fract(f.y));"
   "float p=t(f,vec2(25.8,.5))*.1;"
   "vec3 u=y;"
   "float r=s(i-vec3(2.,5.,0.),vec3(.3));"
   "r=min(r,s(i-vec3(0.,0.,2.),vec3(.3)))-.5-length(sin(y*8.+y.z*15.))*.05;"
   "float x=p;"
   "x=min(x,r*.1);"
   "if(abs(p-x)<.1)"
     "l=1.;"
   "if(abs(r-x)<.1)"
     "l=2.;"
   "return vec2(x*.8,l);"
 "}"
 "vec3 m(vec3 l,vec3 v,vec3 u,vec3 y,vec3 a,float x)"
 "{"
   "return y*(max(0.,dot(v,-u))*.5+max(0.,dot(a,-u))+max(0.,dot(-a,-u))+.1)+pow(max(0.,dot(reflect(v,-u),a))+max(0.,dot(reflect(v,-u),-a)),5.)*x;"
 "}"
 "vec3 m(vec3 l,vec3 a,float y)"
 "{"
   "vec3 u=normal(m,l),f=vec3(1.-step(.001,abs(y-1.)))*.3;"
   "f=f*vec3(cos(v*5.),0.);"
   "return m(l,a,u,f,vec3(0.,0.,1.),1.);"
 "}"
 "vec3 u(vec3 v,vec3 l,float y)"
 "{"
   "vec3 a=normal(x,v);"
   "float u=1.-step(.001,abs(y-1.));"
   "vec3 i=vec3(.3,cos(v.x)*.3,sin(v.y*.5)*.2),p=vec3(abs(sin(f*.5)))*.3,t=mix(i,p,u);"
   "return m(v,l,a,t,vec3(0.,1.,0.),1.);"
 "}"
 "vec3 x(vec3 l,vec3 a,float y)"
 "{"
   "vec3 f=normal(u,l),i=vec3(1.-step(.001,abs(y-1.)))*.3;"
   "i*=vec3(sin(v*2.),1.);"
   "i+=vec3(1,.7,.4)*.3;"
   "return m(l,a,f,i,vec3(0.,0.,1.),.5);"
 "}"
 "vec3 m(vec3 v,vec3 l)"
 "{"
   "vec3 y=vec3(0.);"
   "float f,i=0.;"
   "vec2 a;"
   "vec3 u;"
   "for(int r=0;r<150;r++)"
     "{"
       "u=v+f*l;"
       "a=m(u);"
       "if(a.x<.005||f>50.)"
         "break;"
       "f+=a.x;"
       "i+=1.;"
     "}"
   "float x=dot(l.xy,vec2(0.,1.)),p=pow(abs(x),10.)*4.;"
   "vec3 t=vec3(10.,.2,.2)*(mod(pow(abs(x),1.5),.06)+p);"
   "if(a.x<.01)"
     "y=m(u-.01*l,l,a.y),y=mix(y,p*vec3(1.,.8,.7),min(1.,f/500.));"
   "else"
     " y=t;"
   "y+=pow(abs(i/250.),1.5)*(.7+vec3(1.3,1.3,.6));"
   "return y;"
 "}"
 "vec3 p(vec3 v,vec3 l)"
 "{"
   "vec3 y=vec3(0.);"
   "float f=0.,a=0.,i=0.;"
   "vec2 m;"
   "vec3 p;"
   "for(int r=0;r<80;r++)"
     "{"
       "p=v+f*l;"
       "m=x(p);"
       "if(m.x<.005||f>80.)"
         "break;"
       "f+=m.x;"
       "i+=1.;"
     "}"
   "float t=abs(dot(l.xy,vec2(0.,1.))),s=pow(abs(t),10.)*1.5;"
   "vec3 n=vec3(1.5,1.,.5)*(1.2+cos(l*25.)),r=vec3(2.,.7,.3),d=n*(mod(pow(abs(t),3.),.05)+s);"
   "if(m.x<.1)"
     "y=u(p-.005*l,l,m.y),y=mix(y,s*r,min(1.,f/500.));"
   "else"
     " y=d;"
   "y+=pow(abs(i/80.),3.)*r;"
   "return y;"
 "}"
 "vec3 u(vec3 v,vec3 l)"
 "{"
   "vec3 y=vec3(0.);"
   "float f=0.,a=0.,i=0.;"
   "vec2 m;"
   "vec3 p;"
   "for(int r=0;r<220;r++)"
     "{"
       "p=v+f*l;"
       "m=u(p);"
       "if(m.x<.005||f>100.)"
         "break;"
       "f+=m.x;"
       "i+=1.;"
     "}"
   "float t=dot(l.xy,vec2(0.,1.)),s=pow(abs(t),10.)*1.5;"
   "vec3 n=vec3(s);"
   "float r=m.y;"
   "if(m.x<.1)"
     "y=x(p-.005*l,l,m.y),y=mix(y,s*vec3(1.,.9,.8),min(1.,f/70.));"
   "else"
     " y=n;"
   "y+=pow(abs(i/400.),2.)*(1.+vec3(1.3,0,0.))*(1.-step(1.5,r));"
   "return y;"
 "}"
 "vec3 w(float v)"
 "{"
   "v*=.25;"
   "vec3 y=vec3(cos(v*2.)*4.,-cos(v)*15.,sin(v)*15.);"
   "return y;"
 "}"
 "vec3 h(float v)"
 "{"
   "v*=.25;"
   "vec3 y=vec3(cos(v*2.)*3.,-cos(v)*25.,4.-sin(v)*10.);"
   "return y;"
 "}"
 "vec3 r(float v)"
 "{"
   "v*=.1;"
   "vec3 y=vec3(5.,-cos(v)*100.,5.);"
   "return y;"
 "}"
 "vec3 h(in vec2 v,float y)"
 "{"
   "vec3 l=vec3(0.),f,i,s;"
   "float x=time*(1.+step(53.5,time)+step(a,time));"
   "if(y<1.)"
     "f=normalize(vec3(v,.7)),i=w(x),s=normalize(-i),f=n(s,vec3(1.,1.,0.))*f,l=m(i,f);"
   "if(abs(y-2.)<1.)"
     "f=normalize(vec3(v,1.3)),i=h(x),s=normalize(-i),f=n(s,vec3(1.,1.,0.))*f,l=p(i,f);"
   "if(abs(y-3.)<1.)"
     "v*=t(time*50.),f=normalize(vec3(v,.8)),i=r(x),s=normalize(-i),f=n(s,vec3(1.,1.,0.))*f,l=u(i,f);"
   "return l;"
 "}"
 "vec3 d(vec3 v)"
 "{"
   "vec3 y=v;"
   "v*=1.8;"
   "float l=100.;"
   "for(int f=0;f<5;f++)"
     "v=abs(v*3.)/dot(v,v)-1.,l=min(max(abs(v.x),abs(v.y)),l);"
   "l=pow(max(0.,1.-l),5.);"
   "vec3 f=normalize(y);"
   "float a=time*20.;"
   "f.xy*=mat2(cos(a),-sin(a),sin(a),cos(a));"
   "return abs(f)*l*l*2.+dot(v,v)*.0025;"
 "}"
 "vec3 c(vec2 v)"
 "{"
   "float l,y,a=time*.15-2.,u=2.+sin(a*3.+56.)*1.3;"
   "vec3 m=-vec3(vec2(sin(a*2.))*.15,u),f=normalize(vec3(v,3.)),i=vec3(0),x;"
   "vec2 r;"
   "mat2 p;"
   "l=cos(-a),y=sin(-a);"
   "p=mat2(l,-y,y,l);"
   "m.xz*=p;"
   "f.xz*=p;"
   "l=cos(-a*1.5),y=sin(-a*1.5);"
   "p=mat2(l,-y,y,l);"
   "m.yz*=p;"
   "f.yz*=p;"
   "l=cos(a),y=sin(a);"
   "p=mat2(l,-y,y,l);"
   "m.xy*=p;"
   "f.xy*=p;"
   "float t=.0375;"
   "for(float n=1.;n<40.;n++)"
     "{"
       "float w=n*t;"
       "r=s(m,f,w);"
       "x=m+f*r.x;"
       "i+=d(x)*step(0.,r.x)*smoothstep(0.,u,u-r.x+.4);"
       "x=m+f*r.y;"
       "i+=d(x)*step(0.,r.x)*smoothstep(0.,u,u-r.y+.4);"
     "}"
   "return i*.4;"
 "}"
 "vec3 F(vec3 v)"
 "{"
   "vec3 y=v,l=vec3(0.);"
   "float a=1000.,u=1000.,f,p;"
   "for(float i=0.;i<5.;i++)"
     "{"
       "v=abs(v)/dot(v,v)-.5;"
       "f=min(abs(v.x),abs(v.y));"
       "if(a>f)"
         "a=f,u=length(v),l=normalize(abs(v));"
     "}"
   "float m=pow(max(0.,1.-a),80.);"
   "l=m*abs(normalize(3.+l));"
   "l*=1.5-pow(mod(f-y.y+time*.15-.5,1.),.5);"
   "l+=pow(max(0.,1.-u),7.);"
   "return l*min(1.,time);"
 "}"
 "vec3 F(vec2 v,float y)"
 "{"
   "v*=.6;"
   "float l=time*.05;"
   "vec2 a=vec2(0.,-.25);"
   "vec3 m=-vec3(a,2.5+cos(time*.2+2.)*.5),f=normalize(vec3(v,3.)),i=vec3(0.);"
   "mat2 p=t(y),u=t(time*2.23+3.);"
   "for(float r=20.;r<70.;r++)"
     "{"
       "float x=n(m,f,r*.02);"
       "vec3 s=m+f*x;"
       "s.xz*=p;"
       "s.yz*=u;"
       "vec3 w=F(s)*step(0.,x)*(1.-exp(-.03*(r-20.)));"
       "i=max(i,w);"
     "}"
   "return i*i;"
 "}"
 "vec3 c(vec2 v,vec3 y)"
 "{"
   "vec2 l=v;"
   "l-=vec2(.1,0.);"
   "float a=p(l)-.02;"
   "for(int f=0;f<20;f++)"
     "{"
       "vec2 u=s(vec2(float(f),float(f*2)));"
       "float m=length(l*1000.-u);"
       "if(m<10.1)"
         "a*=m;"
     "}"
   "vec3 f=y*(1.-exp(-4.*pow(abs(a),1.2)));"
   "f=mix(f,vec3(1.,.5,.3),1.-smoothstep(0.,.025,abs(a)));"
   "f+=vec3(2.,.5,.3)*mix(f,vec3(.5),.9-smoothstep(0.,.1,abs(a-y.x/3.5)));"
   "return f*min(1.,time);"
 "}"
 "uint i[]=uint[](lh,li,0u,ln,li,lc,le,0u,lp,le,lo,lp,ll,le,0u,la,lt,0u,lr,le,lv,li,ls,li,lo,ln,999u),D[]=uint[](lt,lh,li,ls,0u,lw,la,ls,0u,ls,lp,la,lc,le,0u,lp,lo,lo,0u,0u,999u),o[]=uint[](lc,0u,lu,0u,la,lg,la,li,ln,0u,la,lt,0u,lf,ll,la,ls,lh,lp,la,lr,lt,ly,0u,l2,l0,l1,l9,999u),b[]=uint[](ls,le,lp,lt,0u,l2,l7,0u,l2,l8,0u,l2,l9,0u,0u,la,lr,lg,le,ln,lt,li,ln,la,999u);"
 "float d(vec2 v,vec4 y)"
 "{"
   "v-=y.xy;"
   "vec2 l=y.zw;"
   "return length(v-l*clamp(dot(v,l)/dot(l,l),0.,1.));"
 "}"
 "float F(vec2 v,uint l,float y)"
 "{"
   "v*=1.+length(v-vec2(.5,.7))*.7;"
   "vec4 f=vec4(.5,-.5,1.,0.);"
   "float a=10.;"
   "a=(l&1u)>0u?min(a,d(v,f.wwwx)):a;"
   "a=(l&2u)>0u?min(a,d(v,f.wxwx)):a;"
   "a=(l&4u)>0u?min(a,d(v,f.wzxw)):a;"
   "a=(l&8u)>0u?min(a,d(v,f.xzxw)):a;"
   "a=(l&16u)>0u?min(a,d(v,f.zxwx)):a;"
   "a=(l&32u)>0u?min(a,d(v,f.zwwx)):a;"
   "a=(l&64u)>0u?min(a,d(v,f.xwxw)):a;"
   "a=(l&128u)>0u?min(a,d(v,f.wwxw)):a;"
   "a=(l&256u)>0u?min(a,d(v,f.wxxw)):a;"
   "a=(l&512u)>0u?min(a,d(v,f.xxxw)):a;"
   "a=(l&1024u)>0u?min(a,d(v,f.xwwx)):a;"
   "a=(l&2048u)>0u?min(a,d(v,f.xxwx)):a;"
   "a=(l&4096u)>0u?min(a,d(v,f.wwxx)):a;"
   "a=(l&8192u)>0u?min(a,d(v,f.xxxx)):a;"
   "a=(l&16384u)>0u?min(a,d(v,f.wzxy)):a;"
   "a=(l&32768u)>0u?min(a,d(v,f.xxxy)):a;"
   "a=(l&65536u)>0u?min(a,d(v,f.wxxy)):a;"
   "a=(l&131072u)>0u?min(a,d(v,f.xwxx)):a;"
   "a=(l&262144u)>0u?min(a,d(v,f.wxxx)):a;"
   "a=(l&524288u)>0u?min(a,d(v,f.xzxy)):a;"
   "return a-y;"
 "}"
 "vec3 F(vec2 v,float l,float y,vec3 x)"
 "{"
   "vec2 a=v;"
   "if(l<1.)"
     "a.x+=1.4;"
   "if(l>1.&&l<2.)"
     "a.x+=1.1;"
   "if(l>2.&&l<3.)"
     "a.x+=1.5;"
   "if(l>3.)"
     "a.x+=1.3;"
   "a.x-=y*5.;"
   "a*=10.;"
   "vec3 f=vec3(0.);"
   "uint u;"
   "for(int m=0;m<30;m++)"
     "{"
       "a.y+=sin(float(m)*2.+time*2.)*.2;"
       "if(l<4.)"
         "u=b[m];"
       "if(l<3.)"
         "u=o[m];"
       "if(l<2.)"
         "u=D[m];"
       "if(l<1.)"
         "u=i[m];"
       "if(u==999u)"
         "break;"
       "f+=pow(max(0.,1.-F(a,u,.1)),13.);"
       "a.x-=1.1;"
     "}"
   "return x+f*min(1.,time*.3);"
 "}"
 "void main()"
 "{"
   "vec2 v=pp;"
   "vec3 l=vec3(0.);"
   "v.x*=1.77;"
   "float y=mod(time*.258333+.25,3.),f=mod((time-84.)*.18,4.),u=smoothstep(.8,1.,fract(f)),i=time*(1.+step(a,time)),x=abs(.5-fract((i-53.)*.258333*4.)),m=step(53.,time)*pow(smoothstep(.3,.5,x),2.-step(a,time)*1.5),p=90.*floor(f)+u*90.;"
   "vec3 t=F(v,p*step(10.,time));"
   "if(time<7.5)"
     "l=c(v,t);"
   "if(time>7.7&&time<45.||time>53.5&&time<83.5)"
     "v=pow(abs(v),1.+vec2(m))*sign(v),l=h(v,y)*(1.+m*3.)+smoothstep(78.,83.5,time)*mod(time,.1)/.1;"
   "if(time>45.5&&time<53.5)"
     "l=c(v);"
   "if(time>84.)"
     "l=F(v,f,u,t);"
   "color=vec4(l,1.);"
 "}";

#endif // FRAGMENT_SHADER_INL_
