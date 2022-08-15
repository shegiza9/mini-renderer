#include "render.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include<iostream>
#include<SDL2/SDL.h>
#include<string>
#include<cmath>
#include<limits>



const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

//Starts up SDL and creates window
bool init();

//Loads media
bool loadMedia();

//Frees media and shuts down SDL
void close();


SDL_Window* gWindow = NULL;
SDL_Renderer* gRenderer = NULL;

bool init() {
	bool success = true;

	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		printf("SDL could not initialize! SDL Error: %s\n", SDL_GetError());
		success = false;
	}

	gWindow = SDL_CreateWindow("my renderer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,SCREEN_WIDTH,SCREEN_HEIGHT,SDL_WINDOW_SHOWN);
	if (gWindow == NULL) {
		printf("Window could not be created! SDL Error: %s\n", SDL_GetError());
		success = false;
	}
	else {
		gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
		if (gRenderer == NULL)
		{
			printf("Renderer could not be created! SDL Error: %s\n", SDL_GetError());
			success = false;
		}
		else {
			SDL_SetRenderDrawColor(gRenderer, 0xff, 0xff, 0xff, 0xff);

		}
	}
	return success;

}

bool loadMedia() {
	bool success = true;
	return success;
}

void close() {
	SDL_DestroyRenderer(gRenderer);
	SDL_DestroyWindow(gWindow);
	gRenderer = NULL;
	gWindow = NULL;

	SDL_Quit();
}



void set_lighttrans(Device_render& device, const Vec4f& eye, const Vec4f& at, const Vec4f& up) {
	device.lighttrans.world = device.transform.world;
	device.lighttrans.view = matrix_set_lookat(eye, at, up);
	device.lighttrans.perspective = matris_set_perspective(PI * 0.7f, 1.333f, 1.0f, 500.0f);
	device.lighttrans.tranform_update();
}


void set_world_rotate(Device_render& device, float theta) {
	Matrix m= matrix_set_rotate( 0, 1, 0, theta);
	device.transform.world = m;
	device.lighttrans.world = device.transform.world;
	device.lighttrans.tranform_update();
	device.transform.tranform_update();
	device.light_pv_multi_object_invv = device.lighttrans.pv_matrix * device.transform.view_inv;
}
void set_my_camera(Device_render &device, const Vec4f& eye,const Vec4f& at,const Vec4f& up) {
	device.transform.set_viewdir(eye, at);
	device.transform.view=matrix_set_lookat(eye, at, up);
	device.transform.tranform_update();
	
}


int main(int argc, char* argv[]) {
	//初始化渲染设备
	/*Mesh* light_plane = new Mesh();
	Triangle triangle;
	triangle.v[0] = vertex_list[0]; triangle.v[1] = vertex_list[1]; triangle.v[2] = vertex_list[2];
	light_plane->triangle_list.push_back(&triangle);
	light_plane->render_state = COLOR;
	light_plane->material = NULL;*/

	std::vector<Mesh*>mesh_list;
	
	Device_render myrender(SCREEN_WIDTH, SCREEN_HEIGHT);
	float camera_x = 5.0f;
	
	const char* filename="Lol_Katarina_Default.obj";
	std::cout << "loading:" << filename << std::endl;
	char* basepath = NULL;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename,
		NULL, true);
	if (!warn.empty()) {
		std::cout << "WARN: " << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << "ERR: " << err << std::endl;
	}

	if (!ret) {
		printf("Failed to load/parse .obj.\n");
		return false;
	}
	std::vector<Material*>my_material;
	for (int j = 0; j < materials.size(); ++j) {
		std::cout << materials[j].diffuse_texname << std::endl;
		std::cout << materials[j].specular_texname << std::endl;
		std::cout << materials[j].normal_texname << std::endl;
		Material* m = new Material();
		m->ambient = {0.2f,0.2f,0.2f ,1.0f};
		m->diffuse = { materials[j].diffuse[0] ,materials[j].diffuse[1] ,materials[j].diffuse[2],1.0f };
		m->specular = { materials[j].specular[0] ,materials[j].specular[1] ,materials[j].specular[2] ,1.0f };
		load_png_image(materials[j].diffuse_texname.c_str(), m->diffuse_texture.texture_map, m->diffuse_texture.width, m->diffuse_texture.height);
		load_png_image(materials[j].specular_texname.c_str(), m->specular_texture.texture_map, m->specular_texture.width, m->specular_texture.height);
		load_png_image(materials[j].normal_texname.c_str(), m->normal_texture.texture_map, m->normal_texture.width, m->normal_texture.height);
		my_material.push_back(m);
	}
	std::cout << materials.size() << std::endl;
	/*std::cout << my_material[0]->diffuse_texture->texture_map.size()<<'\t';
	std::cout << my_material[0]->specular_texture->texture_map.size() << '\t';
	std::cout << my_material[0]->normal_texture->texture_map.size() << '\t';*/
	
	
	float min_x= std::numeric_limits<float>::max(), max_x= std::numeric_limits<float>::min();
	float min_y = min_x, max_y = max_x;
	float min_z = min_x, max_z = max_x;
	std::cout << shapes.size() << std::endl;
	
	for (int p = 0; p < shapes.size();++p) {
		Mesh* mmesh = new Mesh();
		std::vector<Triangle*>triangle_list;
		int indices_size = shapes[p].mesh.indices.size();
		for (int i = 0; i < indices_size-2; i += 3) {
			Triangle* t = new Triangle();
			for (int j = 0; j < 3; ++j) {
				/*int m_id = shape.mesh.material_ids[i + j];
				std::cout << materials[m_id].name << std::endl;*/
				tinyobj::index_t k = shapes[p].mesh.indices[i + j];
				min_x = (attrib.vertices[3 * k.vertex_index] < min_x) ? attrib.vertices[3 * k.vertex_index] : min_x;
				max_x= (attrib.vertices[3 * k.vertex_index] > max_x) ? attrib.vertices[3 * k.vertex_index] : max_x;
				min_y = (attrib.vertices[3 * k.vertex_index + 1] < min_y) ? attrib.vertices[3 * k.vertex_index + 1] : min_y;
				max_y = (attrib.vertices[3 * k.vertex_index + 1] > max_y) ? attrib.vertices[3 * k.vertex_index + 1] : max_y;
				min_z = (attrib.vertices[3 * k.vertex_index + 2] < min_z) ? attrib.vertices[3 * k.vertex_index + 2] : min_z;
				max_z = (attrib.vertices[3 * k.vertex_index + 2] > max_z) ? attrib.vertices[3 * k.vertex_index + 2] : max_z;
				t->set_vertex_pos(j,attrib.vertices[3*k.vertex_index],attrib.vertices[3*k.vertex_index+1],attrib.vertices[3*k.vertex_index+2]);
				t->set_vertex_texcord(j, attrib.texcoords[2*k.texcoord_index],attrib.texcoords[2*k.texcoord_index+1]);
				t->set_vertex_normal(j, attrib.normals[3 * k.normal_index], attrib.normals[3 * k.normal_index + 1], attrib.normals[3 * k.normal_index + 2]);
				t->set_vertex_color(j, 0.58f, 0, 0.827f);
				t->set_vertex_rhw(j, 1.0f);
				
			}
			triangle_list.push_back(t);
		}
		mmesh->triangle_list = triangle_list;
		mmesh->material = my_material[p];
		mmesh->render_state = TEXTURE;
		mmesh->is_illum = false;
		mmesh->stencil = true;
		mesh_list.push_back(mmesh);
	}
	
	float x_cent = (min_x + max_x) * 0.5f;
	float y_cent = (min_y + max_y) * 0.5f;
	float z_cent = (min_z + max_z) * 0.5f;
	Vec4f eye, at, up;
	eye = { x_cent,y_cent,max_z+1.0f,1.0f };
	at = { 0,0,1.0f,1.0f };
	up = { 0,1,0,1 };
	float theta = 1.25;
	Vec4f ll1 = { x_cent-2.0f,max_y,max_z+0.5f,1.0f };
	Vec4f lightcenter = { x_cent,min_y+1.0f,z_cent,1.0f };
	Vec4f lightviewat = ll1 - lightcenter;
	set_lighttrans(myrender, ll1, lightviewat, up);
	Vertex_t vlis[8] = {
		{{ll1.x+0.5f,ll1.y+0.5f,ll1.z+0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x + 0.5f,ll1.y + 0.5f,ll1.z - 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x + 0.5f,ll1.y - 0.5f,ll1.z + 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x + 0.5f,ll1.y - 0.5f,ll1.z - 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x - 0.5f,ll1.y + 0.5f,ll1.z + 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x - 0.5f,ll1.y + 0.5f,ll1.z - 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x - 0.5f,ll1.y - 0.5f,ll1.z + 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
		{{ll1.x - 0.5f,ll1.y - 0.5f,ll1.z - 0.5f,1.0f},{1,0},{1.0f,1.0f,1.0f},{0,0,1,1},1.0f},
	};
	Mesh* lightmesh = new Mesh();
	
	Triangle tri1 = { vlis[6],vlis[2],vlis[0] };
	Triangle tri2 = { vlis[0],vlis[4],vlis[6] };
	Triangle tri3 = { vlis[3],vlis[0],vlis[2] };
	Triangle tri4 = { vlis[3],vlis[1],vlis[0] };
	Triangle tri5 = { vlis[0],vlis[1],vlis[5] };
	Triangle tri6 = { vlis[4],vlis[0],vlis[5] };
	Triangle tri7 = { vlis[6],vlis[7],vlis[2] };
	Triangle tri8 = { vlis[2],vlis[7],vlis[3] };
	Triangle tri9 = { vlis[4],vlis[5],vlis[7] };
	Triangle tri10 = { vlis[6],vlis[4],vlis[7] };
	Triangle tri11 = { vlis[5],vlis[1],vlis[3] };
	Triangle tri12 = { vlis[7],vlis[5],vlis[3] };
	std::vector<Triangle*> light_triangle = { &tri1,&tri2,&tri3,&tri4,
	&tri5,&tri6,&tri7,&tri8,
	&tri9,&tri10,&tri11,&tri12 };
	lightmesh->triangle_list = light_triangle;
	lightmesh->material = NULL;
	lightmesh->is_illum = true;
	lightmesh->stencil = false;
	lightmesh->render_state = COLOR;


	Vertex_t floor_m[4] = {
		{{min_x - 20.0f,min_y - 0.05f,max_z + 20.0f,1.0f},{0,0},{1,1,1},{0,1,0,1},1.0f},
		{{max_x + 20.0f,min_y - 0.05f,max_z + 20.0f,1.0f},{0,0},{1,1,1},{0,1,0,1},1.0f},
		{{min_x - 20.0f,min_y - 0.05f,min_z - 20.0f,1.0f},{0,0},{1,1,1},{0,1,0,1},1.0f},
		{{max_x + 20.0f,min_y - 0.05f,min_z - 20.0f,1.0f},{0,0},{1,1,1},{0,1,0,1},1.0f}
	};
	Mesh* floor_mesh = new Mesh();
	Triangle ftri1 = { floor_m[0],floor_m[1],floor_m[2] };
	Triangle ftri2 = { floor_m[1],floor_m[3],floor_m[2] };
	std::vector<Triangle*> flvec = { &ftri1,&ftri2 };
	floor_mesh->triangle_list = flvec;
	floor_mesh->is_illum = false;
	floor_mesh->stencil = false;
	floor_mesh->render_state = COLOR;
	Material* ground = new Material();
	ground->ambient = { 0.2f,0.2f,0.2f,1.0f };
	ground->diffuse = { 0.9f,0.9f,0.9f,1.0f };
	ground->specular = { 0.4f,0.4f,0.4f,1.0f };
	floor_mesh->material = ground;

	if (!init()) {
		std::cout << "Failed to initialize!\n";
	}
	else {
		
		bool quit = false;
		SDL_Event e;
		while (!quit) {

			while (SDL_PollEvent(&e) != 0) {

				if (e.type == SDL_QUIT) {
					quit = true;
				}
				else if (e.type == SDL_KEYDOWN)
				{
					//Select surfaces based on key press
					switch (e.key.keysym.sym)
					{
					case SDLK_UP:
						camera_x += 0.35f;
						break;

					case SDLK_DOWN:
						camera_x -= 0.35f;
						break;

					case SDLK_LEFT:
						theta -= 0.25;
						break;

					case SDLK_RIGHT:
						theta += 0.25;
						break;

					case SDLK_SPACE:
						break;
					default:
						
						break;
					}
				}
			}
			myrender.device_clear();
			eye.z =max_z+camera_x;

			//更新摄像机矩阵和光摄像机矩阵
			set_my_camera(myrender,eye,at,up);
			set_world_rotate(myrender,theta);
			myrender.set_light(ll1);
			//先得到阴影深度
			myrender.set_shadowbuffer(mesh_list);
			myrender.set_mesh_shadow(floor_mesh);
			//渲染
			myrender.draw_mesh_list(mesh_list,true);
			myrender.draw_mesh(floor_mesh,true);
			myrender.draw_world_light(lightmesh);
			//模板测试
			myrender.do_stencil_test(mesh_list);
			//下面的话坐标系线段的函数可能会报错，因为没对线段3D裁剪
			//myrender.draw_coordinates();
			for (int i = 0; i < SCREEN_HEIGHT; ++i) {
				for (int j = 0; j < SCREEN_WIDTH; ++j) {
					//深度图绘制
					/*float k = myrender.shadowbuffer[i * SCREEN_WIDTH + j];
					k /= myrender.max_deep;
					k = k > 1.0f ? 1.0f : k;
					int color = k * 255;
					CMID(color, 0, 255);
					SDL_SetRenderDrawColor(gRenderer,Uint8(color), Uint8(color), Uint8(color),1);
					SDL_RenderDrawPoint(gRenderer, j, i);*/
					UINT32 color = myrender.framebuffer[i * SCREEN_WIDTH + j];
					SDL_SetRenderDrawColor(gRenderer, (0xff << 16 & color) >> 16, (0xff << 8 & color) >> 8, 0xff & color, (0xff << 24 & color) >> 24);
					SDL_RenderDrawPoint(gRenderer, j, i);
				}
			}
			
			/*for (int i = 0; i < myrender.texture_height; ++i) {
				for (int j = 0; j < myrender.texture_width; ++j) {
					IUINT32 color = myrender.texturebuffer[i * myrender.texture_width + j];
					SDL_SetRenderDrawColor(gRenderer, (0xff << 16 & color) >> 16, (0xff << 8 & color) >> 8, 0xff & color, (0xff << 24 & color) >> 24);
					SDL_RenderDrawPoint(gRenderer, j, i);
				}
			}*/
			SDL_RenderPresent(gRenderer);
		}
	}
	delete lightmesh;
	delete ground;
	delete floor_mesh;
	for (int i = 0; i < my_material.size(); ++i) {
		delete my_material[i];
	}
	for (int i = 0; i < mesh_list.size(); ++i) {
		for (int j = 0; j < mesh_list[i]->triangle_list.size(); ++j) {
			delete mesh_list[i]->triangle_list[j];
		}
		delete mesh_list[i];
	}
	close();

	system("pause");
	return 0;
}