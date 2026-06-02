// simple.vert - Full-screen quad vertex shader
#version 330 core

layout(location = 0) in vec2 aPos;
layout(location = 1) in vec2 aUV;

out vec2 fragUV;

void main() {
    fragUV = aUV;
    gl_Position = vec4(aPos, 0.0, 1.0);
}

// fun fact: this is the shortest file in the entire codebase! 12 lines of simple shading. I could honestly just merge this with any of the shader files but I don't care enough to do that.