#version 430

layout (location=0) out vec4 color;
uniform sampler2D inputTexture;
in vec2 p;

void main()
{
	vec2 coords = 0.5 * p + 0.5;
	//color = texture(inputTexture, coords);

	// Glow sample
	//vec2 coord = vTextureCoord;  
	vec2 coord = coords;  
	float weight[7];  
	// Weight map
	weight[0] = 0.9;  
	weight[1] = 0.227027;  
	weight[2] = 0.1945946;  
	weight[3] = 0.1216216;  
	weight[4] = 0.054054;  
	weight[5] = 0.026216;  
	weight[6] = 0.01216;  
	vec3 result = texture(inputTexture, coord).rgb * weight[0]; 
	vec2 tex_offset = 1.0 / vec2(500.0);  
	for (float i=1.0; i< 7.0; i++) { 
		for (float j=1.0; j< 7.0; j++) { 
			vec3 at = texture(inputTexture, coord + vec2(tex_offset.x * i, tex_offset.y * j)).rgb;
			vec3 bt = texture(inputTexture, coord - vec2(tex_offset.x * i, tex_offset.y * j)).rgb;
			if (at. r + at.b + at.g > 2.2) result += at * weight[int(i)] * weight[int(j)];  
			if (bt. r + bt.b + bt.g > 2.2) result += bt * weight[int(i)] * weight[int(j)];
			   
		} 
	} 

	color = vec4(clamp(result, 0.0, 1.0), 1.0);
}