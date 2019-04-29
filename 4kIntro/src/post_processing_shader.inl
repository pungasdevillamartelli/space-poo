/* File generated with Shader Minifier 1.1.6
 * http://www.ctrl-alt-test.fr
 */
#ifndef POST_PROCESSING_SHADER_INL_
# define POST_PROCESSING_SHADER_INL_

const char *post_processing_shader_glsl =
 "#version 430\n"
 "layout(location=0)out vec4 color;"
 "uniform sampler2D inputTexture;"
 "in vec2 p;"
 "void main()"
 "{"
   "vec2 t=.5*p+.5,x=t;"
   "float v[7];"
   "v[0]=.9;"
   "v[1]=.227027;"
   "v[2]=.194595;"
   "v[3]=.121622;"
   "v[4]=.054054;"
   "v[5]=.026216;"
   "v[6]=.01216;"
   "vec3 i=texture(inputTexture,x).xyz*v[0];"
   "vec2 c=1./vec2(500.);"
   "for(float f=1.;f<7.;f++)"
     "{"
       "for(float z=1.;z<7.;z++)"
         "{"
           "vec3 y=texture(inputTexture,x+vec2(c.x*f,c.y*z)).xyz,r=texture(inputTexture,x-vec2(c.x*f,c.y*z)).xyz;"
           "if(y.x+y.z+y.y>2.2)"
             "i+=y*v[int(f)]*v[int(z)];"
           "if(r.x+r.z+r.y>2.2)"
             "i+=r*v[int(f)]*v[int(z)];"
         "}"
     "}"
   "color=vec4(clamp(i,0.,1.),1.);"
 "}";

#endif // POST_PROCESSING_SHADER_INL_
