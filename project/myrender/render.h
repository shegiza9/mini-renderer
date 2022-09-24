#pragma once
#include "mymath.h"
#include<png.h>
#include<memory>
#include<vector>
#include<cmath>
#include<limits>
#include<tuple>
#include<thread>

int load_png_image(const char* filepath, std::vector<UINT32>& texturebuffer, int& width, int& height);

struct Texture {
    std::vector<UINT32> texture_map;
    int width;
    int height;
};
struct Material {
    Material() {
       /* diffuse_texture->texture_map.resize(512, 0);
        diffuse_texture->width = diffuse_texture->height = 512;
        specular_texture->texture_map.resize(512, 0);
        specular_texture->width = specular_texture->height = 512;
        normal_texture->texture_map.resize(512, 0);
        normal_texture->width = normal_texture->height = 512;*/
    }
    Vec4f ambient;
    Vec4f diffuse;
    Vec4f specular;

    Texture diffuse_texture;
    Texture specular_texture;
    Texture normal_texture;
};
class Mesh {
public:

    std::vector<Triangle*>triangle_list;
    Material* material;
    Render_state render_state;
    bool is_illum;
    bool stencil;
};
//WIREFRAME:线框；TEXTURE：纹理贴图；COLOR：渲染颜色

struct Light {
    Vec4f position;
    Vec4f intensity;

    float constant=1.0f;
    float linear=0.14f;
    float quadratic=0.07f;
};
class Device_render {
public:
    Device_render(int w, int h) {
        width = w;
        height = h;
        transform = TransformMatrix(w, h);
        lighttrans = TransformMatrix(w, h);
        framebuffer.resize(w * h);
        zbuffer.resize(w * h);
        shadowbuffer.resize(w * h);
        stenciltest.resize(w * h);
        background_color = 0xffffe0;
        line_color = 0xffffe0;
        
    };
    void set_light(const Vec4f &pos1) {
        
        l1.position = transform.view * pos1;
        //l2.position = transform.view * pos2;
    }
    //设置像素颜色
    void set_pixel(int x, int y,const UINT32& color) {
        if ((x < width) && (y < height)) {
            framebuffer[y * width + x] = color;
        }
    }
    //读取纹理颜色
    std::tuple<int,int,int,int> texture_read(float u, float v, const Texture& texture,bool is_bilinear_interp) {
        float x = u * (texture.width - 1);
        float y = (1 - v) * (texture.height - 1);
        int x1 = std::floor(x);
        int y1 = std::floor(y);
        x1 = CMID(x1, 0, texture.width - 1);
        y1 = CMID(y1, 0, texture.height - 1);
        if (!is_bilinear_interp) {
            UINT32 color1 = texture.texture_map[y1 * texture.width + x1];
            int a = (color1 >> 24) & 0xff;
            int r = std::floor((color1 >> 16) & 0xff);
            int g = std::floor((color1 >> 8) & 0xff);
            int b = std::floor((color1) & 0xff);
            return std::make_tuple(a,r,g,b);
        }
        //双线性插值
        else {
           
            int x2 = std::ceil(x);
            int y2 = std::ceil(y);
            float t1 = x - x1;
            float t2 = y - y1;

            UINT32 color11 = texture.texture_map[y1 * texture.width + x1];
            UINT32 color12 = texture.texture_map[y2 * texture.width + x1];
            UINT32 color21 = texture.texture_map[y1 * texture.width + x2];
            UINT32 color22 = texture.texture_map[y2 * texture.width + x2];

            int a11 = (color11 >> 24) & 0xff;
            int r11 = (color11 >> 16) & 0xff;
            int g11 = (color11 >> 8) & 0xff;
            int b11 = (color11) & 0xff;

            int a12 = (color12 >> 24) & 0xff;
            int r12 = (color12 >> 16) & 0xff;
            int g12 = (color12 >> 8) & 0xff;
            int b12 = (color12) & 0xff;

            int a21 = (color21 >> 24) & 0xff;
            int r21 = (color21 >> 16) & 0xff;
            int g21 = (color21 >> 8) & 0xff;
            int b21 = (color21) & 0xff;

            int a22 = (color22 >> 24) & 0xff;
            int r22 = (color22 >> 16) & 0xff;
            int g22 = (color22 >> 8) & 0xff;
            int b22 = (color22) & 0xff;

            float one_minus_t1_multi_one_minus_t2 = (1 - t1) * (1 - t2);
            float t1_multi_one_minus_t2 = t1 * (1 - t2);
            float one_minus_t1_multi_t2 = (1 - t1) * t2;
            float t1_multi_t2 = t1 * t2;

            int final_a = one_minus_t1_multi_one_minus_t2 * a11 + t1_multi_one_minus_t2 * a21
                + one_minus_t1_multi_t2 * a12 + t1_multi_t2 * a22;
            int final_r = std::floor(one_minus_t1_multi_one_minus_t2 * r11 + t1_multi_one_minus_t2 * r21
                + one_minus_t1_multi_t2 * r12 + t1_multi_t2 * r22);
            int final_g = std::floor(one_minus_t1_multi_one_minus_t2 * g11 + t1_multi_one_minus_t2 * g21
                + one_minus_t1_multi_t2 * g12 + t1_multi_t2 * g22);
            int final_b = std::floor(one_minus_t1_multi_one_minus_t2 * b11 + t1_multi_one_minus_t2 * b21
                + one_minus_t1_multi_t2 * b12 + t1_multi_t2 * b22);
            final_a = CMID(final_a, 0, 255);
            final_r = CMID(final_r, 0, 255);
            final_g = CMID(final_g, 0, 255);
            final_b = CMID(final_b, 0, 255);

            
            return std::make_tuple(final_a,final_r,final_g,final_b);
        }
        
    }
    //清空buffer
    void device_clear() {
        for (int i = 0; i < height; ++i) {
            int co = (1-(float)i / height) * 255;
            for (int j = 0; j < width; ++j) {
                
                framebuffer[i * width + j] = ((co << 16) | (co << 8) | co);
            }
        }
        //std::fill(framebuffer.begin(), framebuffer.end(), 0xffffe0);
        std::fill(zbuffer.begin(), zbuffer.end(), std::numeric_limits<float>::min());
        std::fill(shadowbuffer.begin(), shadowbuffer.end(), std::numeric_limits<float>::max());
        std::fill(stenciltest.begin(), stenciltest.end(), false);
    }

    //计算重心坐标,p1,p2,p3为逆时针,NDC坐标计算，只用x,y维度
    std::tuple<float, float, float> get_barycentric_coord(const Pos& p1,const Pos& p2,const Pos& p3,const Pos& point) {
        float t = (p1.y - p3.y) * (p2.x - p3.x) + (p2.y - p3.y) * (p3.x - p1.x);
        float b1 = (point.y - p3.y) * (p2.x - p3.x) + (p2.y - p3.y) * (p3.x - point.x);
        float b2 = (point.y - p1.y) * (p3.x - p1.x) + (p3.y - p1.y) * (p1.x - point.x);
        float b3 = (point.y - p2.y) * (p1.x - p2.x) + (p1.y - p2.y) * (p2.x - point.x);
        t = 1.0f / t;

        return std::make_tuple(b1 * t, b2 * t, b3 * t);

    }


    //bresenham画线
    //初始rho=0,每次x+1,rho都+一个斜率m
    //如果rho+m>0.5,则y=y+1,rho=rho+m-1
    //rho+m>0.5->2*rho+2*dy>dx
    void draw_line(int x1, int x2, int y1, int y2, UINT32 c) {
        int dx = x2 - x1;
        int dy = y2 - y1;
        int ux = ((dx > 0) << 1) - 1;//如果大于0则为+1，否则为-1
        int uy = ((dy > 0) << 1) - 1;
        int rho = 0;
        int x = x1, y = y1;
        dx = std::abs(dx);
        dy = std::abs(dy);

        if (dx > dy) {
            for (; x != x2 + ux; x += ux) {
                set_pixel(x, y, c);
                rho += dy;
                if ((rho << 1) >= dx) {
                    rho -= dx;
                    y += uy;
                }
            }
        }
        else {
            for (; y != y2 + uy; y += uy) {
                set_pixel(x, y, c);
                rho += dx;
                if ((rho << 1) >= dy) {
                    rho -= dy;
                    x += ux;
                }
            }
        }
    }

    //绘制扫描线
    void draw_scanline( Mesh* mesh,Scanline& scanline, Vertex_t* world_vertex, Screen_points& my_view,bool is_shadowmap=false,bool outlining=false) {
        int x = scanline.x;
        int y = scanline.y;
        int w = scanline.w;
        int index = y * width;
        for (; w > 0; x++, --w) {
            if (x >= 0 && x < width) {
                if (outlining) {
                    if (stenciltest[index + x] == false) {
                        set_pixel(x, y, 0x00FFFF);
                    }
                    continue;
                }
                float rhw = scanline.v.rhw;//rhw=1/z
                //深度测试
                if (rhw > zbuffer[index + x]) {
                    if (mesh->stencil) stenciltest[index + x] = true;
                    float t_w = 1.0f / rhw;
                    zbuffer[index + x] = rhw;
                    //取屏幕坐标
                    Pos point = { (float)x + 0.5f,(float)y + 0.5f,t_w,1.0f };
                    //在2D坐标中获取重心坐标
                    auto bary_centric_coord = get_barycentric_coord(my_view.p1,my_view.p2,my_view.p3,point);
                    float a = std::get<0>(bary_centric_coord);
                    float b = std::get<1>(bary_centric_coord);
                    float c = std::get<2>(bary_centric_coord);
                    
                    Vertex_t vert = compute_vertex(a, b, c, world_vertex, my_view.p1.w, my_view.p2.w, my_view.p3.w);

                    float visibility = 1.0f;
                    if (is_shadowmap) {
                        //PCSS软阴影
                        Pos vpos = vert.pos;
                        Pos my_pos = this->light_pv_multi_object_invv * vpos;
                        Pos cc;
                        lighttrans.transform_homogenize(cc, my_pos);
                        int shadow_x = cc.x;
                        int shadow_y = cc.y;
                        float my_z = my_pos.w;

                        if (shadow_x>=0&&shadow_x<lighttrans.width&&shadow_y>=0&&shadow_y<lighttrans.height) {
                            float shadow = 0.0;
                            float count = 0.005f;

                            //blocker search
                            float avg_blockdepth;
                            int lightwidth = 3;
                            int kernel = (my_z - 1) / my_z * lightwidth;
                            kernel = std::max(kernel, 0);
                            if (kernel == 0) visibility = 0;
                            else {
                                int size = 0;
                                float depth_sum = 0;
                                for (int wi = -kernel; wi <= kernel; ++wi) {
                                    for (int he = -kernel; he <= kernel; ++he) {
                                        int my_index = (shadow_y + he) * lighttrans.width + (shadow_x + wi);
                                        if (my_index < 0 || my_index >= width * height) continue;
                                        float the_depth = shadowbuffer[my_index];
                                        if (my_z > the_depth) {
                                            size++;
                                            depth_sum += the_depth;
                                        }
                                    }
                                }

                                if (size == 0) visibility = 1.0f;
                                else {
                                    avg_blockdepth = depth_sum / size;
                                    //penumbra estimation
                                    kernel = (my_z - avg_blockdepth) / avg_blockdepth * lightwidth;
                                    if (kernel < 0) kernel = 0;
                                    kernel++;
                                    for (int wi = -kernel; wi <= kernel; ++wi) {
                                        for (int he = -kernel; he <= kernel; ++he) {
                                            int my_index = (shadow_y + he) * lighttrans.width + (shadow_x + wi);
                                            if (my_index < 0 || my_index >= width * height) continue;
                                            float the_depth = shadowbuffer[my_index];
                                            if (my_z > the_depth+0.068f) {
                                                count++;
                                                shadow+=1.0f;
                                            }
                                        }
                                    }
                                    shadow = shadow / count;
                                    visibility = 1.0f - shadow;
                                }
                            }
                           
                            /*for (int wi = -2; wi < 3; ++wi) {
                                for (int he = -2; he < 3; ++he) {
                                    int my_index = (shadow_y + he) * lighttrans.width + (shadow_x + wi);
                                    if (my_index<0 || my_index>=width * height) continue;
                                    if (my_z > shadowbuffer[my_index] + 0.068f) {
                                        shadow += 1.0f;
                                        

                                    }
                                    count++;
                                }
                            }*/
                            
                        }
                        
                        
                    }
                    

                    Vec4f la, ld, ls;
                    
                    if (mesh->render_state == COLOR) {
                        if (!mesh->is_illum) {
                            set_frag_shader(mesh, vert, la, ld, ls, NULL, false);
                            Vec4f kk = la + ld * visibility + ls * visibility;
                            float r = vert.color.r * kk.x;
                            float g = vert.color.g * kk.y;
                            float b = vert.color.b * kk.z;
                            int R, G, B;
                            
                
                            R = (int)(r * 255.0f );
                            G = (int)(g * 255.0f );
                            B = (int)(b * 255.0f );
              
                            
                            R = CMID(R, 0, 255);
                            G = CMID(G, 0, 255);
                            B = CMID(B, 0, 255);
                            set_pixel(x, y, (R << 16) | (G << 8) | (B));
                        }
                        else {
                            float r = vert.color.r;
                            float g = vert.color.g;
                            float b = vert.color.b;
                            int R, G, B;
                            R = (int)(r * 255.0f * visibility);
                            G = (int)(g * 255.0f * visibility);
                            B = (int)(b * 255.0f * visibility);
                            R = CMID(R, 0, 255);
                            G = CMID(G, 0, 255);
                            B = CMID(B, 0, 255);
                            set_pixel(x, y, (R << 16) | (G << 8) | (B));
                        }
                    }
                    if (mesh->render_state == TEXTURE) {
                        float u = vert.tc.u;
                        float v = vert.tc.v;
                        std::tuple<int, int, int, int> argb_dif=texture_read(u, v,mesh->material->diffuse_texture,false);
                        std::tuple<int, int, int, int> argb_spec = texture_read(u, v, mesh->material->specular_texture, false);
                        std::tuple<int, int, int, int> argb_norm = texture_read(u, v, mesh->material->normal_texture, false);
                        int R1 = std::get<1>(argb_dif);
                        int G1 = std::get<2>(argb_dif);
                        int B1 = std::get<3>(argb_dif);

                        int R2 = std::get<1>(argb_spec);
                        int G2 = std::get<2>(argb_spec);
                        int B2 = std::get<3>(argb_spec);

                        int R3 = std::get<1>(argb_norm);
                        int G3 = std::get<2>(argb_norm);
                        int B3 = std::get<3>(argb_norm);

                        float nx = (float)R3 / 255;
                        float ny = (float)G3 / 255;
                        float nz = (float)B3 / 255;

                        Vec4f m_normal = { nx*2.0f-1.0f,ny*2.0f-1.0f,nz*2.0f-1.0f,1.0f };
                        set_frag_shader(mesh, vert, la, ld, ls, &m_normal, false);
                        
                        int R, G, B;
                        
                        ld = ld * visibility;
                        ls = ls * visibility;
                        R = (R1 * (la.x + ld.x) + R2 * ls.x) ;
                        G = (G1 * (la.y + ld.y) + G2 * ls.y) ;
                        B = (B1 * (la.z + ld.z) + B2 * ls.z) ;
                        
                        
                        

                        /*int R = std::pow((R1 * (la.x + ld.x) + R2 * ls.x)/255,1/gamma)*255;
                        int G = std::pow((G1 * (la.y + ld.y) + G2 * ls.y)/255,1/gamma)*255;
                        int B = std::pow((B1 * (la.z + ld.z) + B2 * ls.z)/255,1/gamma)*255;*/

                        /*int R = R1 * (la.x + ld.x + ls.x);
                        int G = G1 * (la.y + ld.y + ls.y);
                        int B = B1 * (la.z + ld.z + ls.z);*/

                        R = CMID(R, 0, 255);
                        G = CMID(G, 0, 255);
                        B = CMID(B, 0, 255);
                        set_pixel(x, y, (R << 16) | (G << 8) | (B));
                    }
                }
            }
            vertex_add(scanline.v, scanline.step);
            if (x >= width) break;
        }
    }

    void render_shadow_scanline(Scanline& scanline) {
        int x = scanline.x;
        int y = scanline.y;
        int w = scanline.w;
        int index = y * lighttrans.width;
        for (; w > 0; ++x, --w) {
            if (x >= 0 && x < lighttrans.width) {
                float rhw = 1.0f/scanline.v.rhw;//rhw=1/z;
                if (rhw < shadowbuffer[index + x]) {
                    if (rhw < 1.0f) rhw = 1.0f;
                    shadowbuffer[index + x] = rhw;
                    if (rhw > max_deep) max_deep = rhw;
                }
            }
            vertex_add(scanline.v, scanline.step);
            if (x >= width) break;
        }
    }
    //渲染平底三角形
    void render_trapezoid( Mesh* mesh,Trapezoid& trape,Vertex_t* world_vertex,Screen_points &my_view,bool is_shadowmap,bool outlining) {
        int top = (int)(trape.top + 0.5f);
        int bottom = (int)(trape.bottom + 0.5f);
        Scanline scan;
        for (int j = top; j < bottom; ++j) {
            if (j >= 0 && j < height) {
                trapezoid_edge_interp(trape, (float)j + 0.5f);
                trapezoid_init_scan_line(trape, scan, j);
                draw_scanline(mesh,scan,world_vertex,my_view,is_shadowmap,outlining);
            }
            if (j >= height) break;
        }
    }
    void render_shadow_trapezoid(Trapezoid& trape) {
        int top = (int)(trape.top + 0.5f);
        int bottom = (int)(trape.bottom + 0.5f);
        Scanline scan;
        for (int j = top; j < bottom; ++j) {
            if (j >= 0 && j < height) {
                trapezoid_edge_interp(trape, (float)j + 0.5f);
                trapezoid_init_scan_line(trape, scan, j);
                render_shadow_scanline(scan);
            }
            if (j >= height) break;
        }
    }
    
    void render_shadow_triangle(Triangle* tri) {
        Vertex_t t1 = tri->v[0], t2 = tri->v[1], t3 = tri->v[2];
        Pos p1, p2, p3, c1, c2, c3;
        p1 = lighttrans.transform * t1.pos;
        p2 = lighttrans.transform * t2.pos;
        p3 = lighttrans.transform * t3.pos;

        lighttrans.transform_homogenize(c1, p1);
        lighttrans.transform_homogenize(c2, p2);
        lighttrans.transform_homogenize(c3, p3);

        t1.pos = c1;
        t2.pos = c2;
        t3.pos = c3;

        t1.pos.w = p1.w;
        t2.pos.w = p2.w;
        t3.pos.w = p3.w;

        vertex_init_rhw(t1);
        vertex_init_rhw(t2);
        vertex_init_rhw(t3);

        Trapezoid trap[2];
        int len = trapezoid_init_triangle(trap, &t1, &t2, &t3);
        if (len >= 1) render_shadow_trapezoid(trap[0]);
        if (len >= 2) render_shadow_trapezoid(trap[1]);
    }
    //渲染任意三角形
    void render_primitive( Mesh* mesh,Triangle* tri,bool is_shadowmap) {
        Vertex_t t1 = tri->v[0], t2 = tri->v[1], t3 = tri->v[2];
        //先转换为视角坐标
        //Pos p1, p2, p3, c1, c2, c3;

        
        t1.pos = transform.vm_matrix * t1.pos;
        t2.pos = transform.vm_matrix * t2.pos;
        t3.pos = transform.vm_matrix* t3.pos;

        t1.tangent = transform.vm_matrix * t1.tangent;
        t2.tangent = transform.vm_matrix * t2.tangent;
        t3.tangent = transform.vm_matrix * t3.tangent;

        //顶点向量也跟着转换一下
        Matrix n = transform.vm_inv_tran;
        t1.normal = n * t1.normal;
        t2.normal = n * t2.normal;
        t3.normal = n * t3.normal;
        t1.normal.normalized();
        t2.normal.normalized();
        t3.normal.normalized();
        //背面剔除
        

        ////gouraud着色
        //set_gouraud_shader(mesh, t1);
        //set_gouraud_shader(mesh, t2);
        //set_gouraud_shader(mesh, t3);
        // 
        if (is_back(t1.pos,t2.pos,t3.pos)) return;
        
        vertex_culling(t1, t2, t3, mesh,is_shadowmap,false);
        //render_triangle(t1, t2, t3, mesh);


    }
    
    void render_enlarged_primitive(Mesh* mesh, Triangle* tri) {
        Vertex_t t1 = tri->v[0], t2 = tri->v[1], t3 = tri->v[2];
        //先转换为视角坐标
        t1.pos = t1.pos + t1.normal * 0.05f;
        t2.pos = t2.pos + t2.normal * 0.05f;
        t3.pos = t3.pos + t3.normal * 0.05f;

        t1.pos = transform.vm_matrix * t1.pos;
        t2.pos = transform.vm_matrix * t2.pos;
        t3.pos = transform.vm_matrix * t3.pos;
        if (is_back(t1.pos, t2.pos, t3.pos)) return;

        vertex_culling(t1, t2, t3, mesh,false,true);
    }
    //t1,t2,t3为视角坐标顶点，mesh为Mesh指针，传递material参数
    void render_triangle(Vertex_t &tt1,Vertex_t &tt2,Vertex_t &tt3,Mesh* mesh,bool is_shadowmap,bool outlining) {
        
        Pos p1, p2, p3, c1, c2, c3;
        Vertex_t t1 = tt1, t2 = tt2, t3 = tt3;
        //储存视角坐标,法线，颜色,为phong着色做准备
        Vertex_t world[3] = { t1,t2,t3 };
        Vec4f tangent, bitangent;
        

        
        //储存pvm后的坐标
        c1 = transform.perspective * t1.pos;
        c2 = transform.perspective * t2.pos;
        c3 = transform.perspective * t3.pos;


        //屏幕裁剪
        int check1 = transform.transform_check_cvv(c1);
        int check2 = transform.transform_check_cvv(c2);
        int check3 = transform.transform_check_cvv(c3);

        //裁减图元
        // z<-w:check|=1;z>w:check|=2;
        //x<-w:check|=4;x>w:check|=8;
        //y<-w;check|=16;y>w:check|=32
        //if ((check1 & 1) && (check2 & 1) && (check3 & 1)) return;
        //if ((check1 & 2) && (check2 & 2) && (check3 & 2)) return;

        if ((check1 & 4) && (check2 & 4) && (check3 & 4)) return;
        if ((check1 & 8) && (check2 & 8) && (check3 & 8)) return;

        if ((check1 & 16) && (check2 & 16) && (check3 & 16)) return;
        if ((check1 & 32) && (check2 & 32) && (check3 & 32)) return;





        //变为屏幕坐标
        transform.transform_homogenize(p1, c1);
        transform.transform_homogenize(p2, c2);
        transform.transform_homogenize(p3, c3);
        Screen_points my_view;
        my_view.p1 = p1;
        my_view.p2 = p2;
        my_view.p3 = p3;
        my_view.p1.w = 1.0f / c1.w;
        my_view.p2.w = 1.0f / c2.w;
        my_view.p3.w = 1.0f / c3.w;

        if (mesh->render_state == TEXTURE || mesh->render_state == COLOR) {
            Trapezoid trap[2];

            t1.pos = p1;
            t2.pos = p2;
            t3.pos = p3;
            t1.pos.w = c1.w;
            t2.pos.w = c2.w;
            t3.pos.w = c3.w;


            vertex_init_rhw(t1);
            vertex_init_rhw(t2);
            vertex_init_rhw(t3);

            int len = trapezoid_init_triangle(trap, &t1, &t2, &t3);
            if (len >= 1) render_trapezoid(mesh, trap[0], world, my_view,is_shadowmap,outlining);
            if (len >= 2) render_trapezoid(mesh, trap[1], world, my_view,is_shadowmap,outlining);
        }
        if (mesh->render_state == WIREFRAME) {
            draw_line(int(p1.x), int(p2.x), int(p1.y), int(p2.y), line_color);
            draw_line(int(p1.x), int(p3.x), int(p1.y), int(p3.y), line_color);
            draw_line(int(p2.x), int(p3.x), int(p2.y), int(p3.y), line_color);
        }
        

    }
    

    Vertex_t compute_vertex(float a, float b, float c, Vertex_t* world_vertex, float w1, float w2, float w3) {
        float a_mul_invw1 = a * w1;
        float b_mul_invw2 = b * w2;
        float c_mul_invw3 = c * w3;
        float interp_invw = a_mul_invw1 + b_mul_invw2 + c_mul_invw3;
        float inv_interp_invw = 1.0f / interp_invw;
        Vertex_t vert;
        vert.pos = world_vertex[0].pos * a_mul_invw1 + world_vertex[1].pos * b_mul_invw2 + world_vertex[2].pos * c_mul_invw3;
        vert.pos = vert.pos * inv_interp_invw;

        vert.normal = world_vertex[0].normal * a_mul_invw1 + world_vertex[1].normal * b_mul_invw2 + world_vertex[2].normal * c_mul_invw3;
        vert.normal = vert.normal * inv_interp_invw;
        vert.normal.normalized();

        vert.tangent = world_vertex[0].tangent * a_mul_invw1 + world_vertex[1].tangent * b_mul_invw2 + world_vertex[2].tangent * c_mul_invw3;
        vert.tangent = vert.tangent * inv_interp_invw;
        vert.tangent.normalized();

        vert.color = world_vertex[0].color * a_mul_invw1 + world_vertex[1].color * b_mul_invw2 + world_vertex[2].color * c_mul_invw3;
        vert.color = vert.color * inv_interp_invw;

        vert.tc.u = world_vertex[0].tc.u * a_mul_invw1 + world_vertex[1].tc.u * b_mul_invw2 + world_vertex[2].tc.u * c_mul_invw3;
        vert.tc.u *= inv_interp_invw;
        vert.tc.v = world_vertex[0].tc.v * a_mul_invw1 + world_vertex[1].tc.v * b_mul_invw2 + world_vertex[2].tc.v * c_mul_invw3;
        vert.tc.v *= inv_interp_invw;
        return vert;
    }

    
    void draw_mesh_list(std::vector<Mesh*>& mesh_list,bool is_shadowmap) {
        for (auto mesh : mesh_list) {
            draw_mesh(mesh,is_shadowmap);
        }
    }
    void draw_mesh(Mesh* mesh,bool is_shadowmap) {
        /*for (int i = 0; i < mesh->triangle_list.size(); ++i) {
            render_primitive(mesh, mesh->triangle_list[i],is_shadowmap);
        }*/
        int len = mesh->triangle_list.size();

        std::thread thread1(&Device_render::draw_mesh_thread, this, mesh, is_shadowmap, 0, len / 4);
        std::thread thread2(&Device_render::draw_mesh_thread, this, mesh, is_shadowmap, len / 4, len / 2);
        std::thread thread3(&Device_render::draw_mesh_thread, this, mesh, is_shadowmap, len / 2, len * 3 / 4);
        std::thread thread4(&Device_render::draw_mesh_thread, this, mesh, is_shadowmap, len * 3 / 4, len);
        thread1.join();
        thread2.join();
        thread3.join();
        thread4.join();

       
        
    }
    void draw_mesh_thread(Mesh* mesh, bool is_shadowmap, int begin, int end) {
        for (int i = begin; i < end; ++i) {
            render_primitive(mesh, mesh->triangle_list[i], is_shadowmap);
        }
    }

    void set_shadowbuffer(std::vector<Mesh*>& mesh_list) {
        for (auto mesh : mesh_list) {
            set_mesh_shadow(mesh);
        }
    }
    void set_mesh_shadow(Mesh* mesh) {

        /*for (int i = 0; i < mesh->triangle_list.size(); ++i) {
            render_shadow_triangle(mesh->triangle_list[i]);
        }*/
        
        int len = mesh->triangle_list.size();

        std::thread thread1(&Device_render::set_shadow_thread, this, mesh, 0, len / 4);
        std::thread thread2(&Device_render::set_shadow_thread, this, mesh, len / 4, len / 2);
        std::thread thread3(&Device_render::set_shadow_thread, this, mesh, len / 2, len / 4 * 3);
        std::thread thread4(&Device_render::set_shadow_thread, this, mesh, 0, len);
        thread1.join();
        thread2.join();
        thread3.join();
        thread4.join();

        /*
        else {
            std::thread thread1(&Device_render::set_shadow_thread, this, mesh, 0, len);
            thread1.join();
        }
        */
    }
    void set_shadow_thread(Mesh* mesh, int begin, int end) {
        for (int i = begin; i < end; ++i) {
            render_shadow_triangle(mesh->triangle_list[i]);
        }
    }
    void do_stencil_test(std::vector<Mesh*>& mesh_list) {
        for (auto mesh : mesh_list) {
            for (int i = 0; i < mesh->triangle_list.size(); ++i) {
                render_enlarged_primitive(mesh, mesh->triangle_list[i]);
            }
        }
    }
    void draw_world_light(Mesh* mesh) {
        for (int i = 0; i < mesh->triangle_list.size(); ++i) {
            Vertex_t p1 = mesh->triangle_list[i]->v[0], p2 = mesh->triangle_list[i]->v[1], p3 = mesh->triangle_list[i]->v[2];
            p1.pos = transform.view * p1.pos;
            p2.pos = transform.view * p2.pos;
            p3.pos = transform.view * p3.pos;
            render_triangle(p1, p2, p3, mesh,false,false);
        }
    }
    void draw_coordinates() {
        draw_light_point(l1.position);
        draw_axis_xyz({ 2.0f,0,0,1.0f }, { -2.0f,0,0,1.0f }, line_color);
        draw_axis_xyz({ 0,1.0f,0,1.0f }, { 0,-1.0f,0,1.0f }, 0xB22222);
        draw_axis_xyz({ 0,0,1.0f,1.0f }, { 0,0,-1.0f,1.0f }, 0x00FF7F);
    }
    void draw_light_point(const Pos &p) {
        Vec4f light1 = transform.perspective * p;
        

        Vec4f ml1;
        transform.transform_homogenize(ml1, light1);

        set_pixel(ml1.x, ml1.y, 0xffffff);
        //set_pixel(ml2.x, ml2.y, 0xffffff);
    }

    void draw_axis_xyz(const Vec4f& x11,const Vec4f& x22,UINT32 cc) {

        Vec4f mx11 = transform.transform * x11;
        Vec4f mx22 = transform.transform * x22;
        Vec4f hx11, hx22;
        transform.transform_homogenize(hx11, mx11);
        transform.transform_homogenize(hx22, mx22);
        draw_line(int(hx11.x), int(hx22.x), int(hx11.y), int(hx22.y), cc);
    }

    //背面剔除,true:背面，false:正面,p1,p2,p3为逆时针
    bool is_back(Pos &p1,Pos &p2,Pos &p3) {
        Vec4f edge1 = p2 - p1;
        Vec4f edge2 = p3 - p1;

        Vec4f normal = edge1.crossProduct(edge2);
        Vec4f eye(0, 0, 0, 1.0f);
        Vec4f viewdir = p1 - eye;
        

        if (normal.dotProduct(viewdir) > 0.01f) return true;
        else return false;
    }
    //只对近z面裁剪
    void vertex_culling(Vertex_t &t1,Vertex_t &t2,Vertex_t &t3,Mesh* mesh,bool is_shadowmap,bool outlining) {
        int vertex_codes[3]={0,0,0};
        //0:在nearz里，1:在nearz外
        int num_vertex_in = 3;
        float nearz = -transform.nz;
        if (t1.pos.z > nearz) {
            vertex_codes[0] = 1;
            num_vertex_in--;
        }
        if (t2.pos.z > nearz) {
            vertex_codes[1] = 1;
            num_vertex_in--;
        }
        if (t3.pos.z > nearz) {
            vertex_codes[2] = 1;
            num_vertex_in--;
        }
        
        Vertex_t temp;
        Vertex_t vi,vj;
        
        //全部在里面，不用裁剪
        if (num_vertex_in == 3) {
            render_triangle(t1, t2, t3, mesh,is_shadowmap,outlining);
            return;
        }
        //全部在外面，直接舍去
        else if (num_vertex_in == 0) {
            //std::cout << "out" <<std::endl;
            return;
        }
        //只有一个顶点在里面
        else if (num_vertex_in == 1) {
            //交换顺序，使得逆时针情况下p1在里面
            if (vertex_codes[0] == 0) {

            }
            else if (vertex_codes[1] == 0) {
                temp = t1;
                t1 = t2;
                t2 = t3;
                t3 = temp;
            }
            else {
                temp = t1;
                t1 = t3;
                t3 = t2;
                t2 = temp;
            }

            float tt1 = (nearz - t1.pos.z) / (t2.pos.z - t1.pos.z);
            float tt2 = (nearz - t1.pos.z) / (t3.pos.z - t1.pos.z);

            float xi = interp(t1.pos.x, t2.pos.x, tt1);
            float yi = interp(t1.pos.y, t2.pos.y, tt1);

            t2.pos.x = xi;
            t2.pos.y = yi;
            t2.pos.z = nearz;

            float xj = interp(t1.pos.x, t3.pos.x, tt2);
            float yj = interp(t1.pos.y, t3.pos.y, tt2);

            t3.pos.x = xj;
            t3.pos.y = yj;
            t3.pos.z = nearz;

            t2.color = t1.color * (1 - tt1) + t2.color * tt1;
            t2.normal = t1.normal.vec_interpolate(t2.normal, tt1);
            t2.normal.normalized();
            t2.tangent = t1.tangent.vec_interpolate(t2.tangent, tt1);
            t2.tangent.normalized();

            t2.tc.u = t1.tc.u * (1 - tt1) + t2.tc.u * tt1;
            t2.tc.v = t1.tc.v * (1 - tt1) + t2.tc.v * tt1;

            t3.color = t1.color * (1 - tt2) + t3.color * tt2;
            t3.normal = t1.normal.vec_interpolate(t3.normal, tt2);
            t3.normal.normalized();
            t3.tangent = t1.tangent.vec_interpolate(t3.tangent, tt1);
            t3.tangent.normalized();

            t3.tc.u = t1.tc.u * (1 - tt2) + t3.tc.u * tt2;
            t3.tc.v = t1.tc.v * (1 - tt2) + t3.tc.v * tt2;

            render_triangle(t1, t2, t3, mesh,is_shadowmap,outlining);

        }
        //有两个顶点在视锥体内
        else {
            if (vertex_codes[0] == 1) {

            }
            else if (vertex_codes[1] == 1) {
                temp = t1;
                t1 = t2;
                t2 = t3;
                t3 = temp;
            }
            else {
                temp = t1;
                t1 = t3;
                t3 = t2;
                t2 = temp;
            }

            float tt1 = (nearz - t1.pos.z) / (t2.pos.z - t1.pos.z);
            float tt2 = (nearz - t1.pos.z) / (t3.pos.z - t1.pos.z);

            float xi = t1.pos.x * (1 - tt1) + t2.pos.x * tt1;
            float yi = t1.pos.y * (1 - tt1) + t2.pos.y * tt1;
            
            
            
            vi.pos = { xi,yi,nearz,1.0f };
            vi.color = t1.color * (1 - tt1) + t2.color * tt1;
            vi.normal = t1.normal * (1 - tt1) + t2.normal * tt1;
            vi.normal.normalized();

            vi.tangent = t1.tangent * (1 - tt1) + t2.tangent * tt1;
            vi.tangent.normalized();

            vi.tc.u = t1.tc.u * (1 - tt1) + t2.tc.u * tt1;
            vi.tc.v = t1.tc.v * (1 - tt1) + t2.tc.v * tt1;
            
            
            float xj = t1.pos.x * (1 - tt2) + t3.pos.x * tt2;
            float yj = t1.pos.y * (1 - tt2) + t3.pos.y * tt2;

            vj.pos = { xj,yj,nearz,1.0f };
            vj.color = t1.color * (1 - tt2) + t3.color * tt2;
            vj.normal = t1.normal * (1 - tt2) + t3.normal * tt2;
            vj.normal.normalized();

            vj.tangent = t1.tangent * (1 - tt2) + t3.tangent * tt2;
            vj.tangent.normalized();

            vj.tc.u = t1.tc.u * (1 - tt2) + t3.tc.u * tt2;
            vj.tc.v = t1.tc.v * (1 - tt2) + t3.tc.v * tt2;

            render_triangle(vi, t2, t3, mesh,is_shadowmap,outlining);
            render_triangle(vi, t3, vj, mesh,is_shadowmap,outlining);
            //render_triangle(p1, p2, p3, mesh);
        }
    }

    
    //相机坐标下进行
    void set_frag_shader(Mesh* mesh, Vertex_t& v,Vec4f &kia,Vec4f &kid,Vec4f &kis,Vec4f* m_normal,bool is_normal) {

        /*Light l1 = { {-1.0f,5.0f,1.0f,1.0f},{1.0f,1.0f,1.0f,1.0f} };
        Light l2 = { {-3.0f,0.0f,1.5f,1.0f},{1.0f,1.0f,1.0f,1.0f} };*/
        /*Light l1 = { {1.0,0,-5.0f,1.0f},{1.0f,1.0f,1.0f,1.0f} };*/
        
        std::vector<Light>lights={l1};
        
        //环境光
        Vec4f ambient_intensity = {0.2f,0.2f,0.2f,1.0f};
        Vec4f La = mesh->material->ambient ;

        //diffuse 
        Vec4f Ld = {0,0,0,1.0f};
        Vec4f Ls = { 0,0,0,1.0f };

        Vec4f eye = { 0,0,0,1.0f };
        Vec4f norm=v.normal;
        if (is_normal) {
            Vec4f t = v.tangent;
            t = t - norm * (norm.dotProduct(t));
            Vec4f bi = norm.crossProduct(t);
            Matrix trans = { t.x,bi.x,norm.x,0.0f,
            t.y,bi.y,norm.y,0.0f,
            t.z,bi.z,norm.z,0.0f,
            0.0f,0.0f,0.0f,1.0f };
            norm = trans * (*m_normal);
        }
        for (auto light : lights) {
            Vec4f pos_to_eye=eye-v.pos;
            pos_to_eye.normalized();

            Vec4f v_to_light = light.position - v.pos;
            float r_1 = v_to_light.length();
            float r_2 = v_to_light.dotProduct(v_to_light);
            v_to_light.normalized();

            Vec4f add_h = (pos_to_eye + v_to_light).normalize();
            float attenuation = 1.0f / (light.constant + light.linear * r_1 + light.quadratic * r_2);
            Vec4f L_d = light.intensity * mesh->material->diffuse;
            L_d = L_d * attenuation * std::max(0.0f, v_to_light.dotProduct(norm));
            Ld += L_d;

            Vec4f L_s = light.intensity * mesh->material->specular;
            L_s = L_s * attenuation * std::pow(std::max(0.0f, norm.dotProduct(add_h)), 10);
            Ls += L_s;
        }


        kia = La;
        kid = Ld;
        kis = Ls;
         
    }
    TransformMatrix transform;
    TransformMatrix lighttrans;
    Matrix light_pv_multi_object_invv;
    int width;
    int height;
    std::vector<UINT32> framebuffer;
    std::vector<float> zbuffer;
    std::vector<float> shadowbuffer;
    std::vector<bool> stenciltest;
    UINT32 background_color;
    UINT32 line_color;
    float max_deep= std::numeric_limits<float>::min();
    Light l1 = { {-1.0f,8.0f,1.0f,1.0f},{1.0f,1.0f,1.0f,1.0f} };
    Light l2 = { {-3.0f,0.0f,1.5f,1.0f},{1.0f,1.0f,1.0f,1.0f} };
};

#define PNG_BYTES_TO_CHECK 8
int load_png_image(const char* filepath, std::vector<UINT32>& texturebuffer, int& width, int& height)
{
    FILE* fp;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep* row_pointers;
    char buf[PNG_BYTES_TO_CHECK];
    int w, h, x, y, temp, color_type;

    errno_t err = fopen_s(&fp, filepath, "rb");
    if (err != 0) {
        return 1; /* 返回值 */
    }

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    info_ptr = png_create_info_struct(png_ptr);

    setjmp(png_jmpbuf(png_ptr));
    /* 读取PNG_BYTES_TO_CHECK个字节的数据 */
    temp = (int)fread(buf, 1, PNG_BYTES_TO_CHECK, fp);
    /* 若读到的数据并没有PNG_BYTES_TO_CHECK个字节 */
    if (temp < PNG_BYTES_TO_CHECK) {
        fclose(fp);
        png_destroy_read_struct(&png_ptr, &info_ptr, 0);
        return 2;/* 返回值 */
    }
    /* 检测数据是否为PNG的签名 */
    temp = png_sig_cmp((png_bytep)buf, (png_size_t)0, PNG_BYTES_TO_CHECK);
    /* 如果不是PNG的签名，则说明该文件不是PNG文件 */
    if (temp != 0) {
        fclose(fp);
        png_destroy_read_struct(&png_ptr, &info_ptr, 0);
        return 3;/* 返回值 */
    }

    /* 复位文件指针 */
    rewind(fp);
    /* 开始读文件 */
    png_init_io(png_ptr, fp);
    /* 读取PNG图片信息 */
    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_EXPAND, 0);
    /* 获取图像的色彩类型 */
    color_type = png_get_color_type(png_ptr, info_ptr);
    /* 获取图像的宽高 */
    w = png_get_image_width(png_ptr, info_ptr);
    h = png_get_image_height(png_ptr, info_ptr);

    texturebuffer.resize(w * h);
    /* 获取图像的所有行像素数据，row_pointers里边就是rgba数据 */
    row_pointers = png_get_rows(png_ptr, info_ptr);
    /* 根据不同的色彩类型进行相应处理 */
    switch (color_type) {
    case PNG_COLOR_TYPE_RGB_ALPHA:
    {
        for (y = 0; y < h; ++y) {
            for (x = 0; x < w; ++x) {
                texturebuffer[y * w + x] = 0;
                /* 以下是RGBA数据，需要自己补充代码，保存RGBA数据 */
                texturebuffer[y * w + x] |= row_pointers[y][4 * x + 0] << 16; // red
                texturebuffer[y * w + x] |= row_pointers[y][4 * x + 1] << 8; // green
                texturebuffer[y * w + x] |= row_pointers[y][4 * x + 2]; // blue
                texturebuffer[y * w + x] |= row_pointers[y][4 * x + 3] << 24; // alpha
            }
        }
        std::cout << "rgba" << std::endl;
    }
    break;

    case PNG_COLOR_TYPE_RGB:
    {
        for (y = 0; y < h; ++y) {
            for (x = 0; x < w; ++x) {
                texturebuffer[y * w + x] = 0xff000000;
                texturebuffer[y * w + x] |= row_pointers[y][3 * x + 0] << 16; // red
                texturebuffer[y * w + x] |= row_pointers[y][3 * x + 1] << 8; // green
                texturebuffer[y * w + x] |= row_pointers[y][3 * x + 2]; // blue
            }
        }
        std::cout << "rgb" << std::endl;
    }
    break;
    case PNG_COLOR_TYPE_GRAY:
    {
        for (y = 0; y < h; ++y) {
            for (x = 0; x < w; ++x) {
                texturebuffer[y * w + x] = 0xff000000;
                texturebuffer[y * w + x] |= row_pointers[y][x] << 16; 
                texturebuffer[y * w + x] |= row_pointers[y][x] << 8; 
                texturebuffer[y * w + x] |= row_pointers[y][x]; 
            }
        }
        std::cout << "gray" << std::endl;
    }
    break;
    default:
        fclose(fp);
        png_destroy_read_struct(&png_ptr, &info_ptr, 0);
        std::cout << "other" << std::endl;
        return 4/* 返回值 */;
    }
    png_destroy_read_struct(&png_ptr, &info_ptr, 0);
    std::cout << w << std::endl;
    std::cout << h << std::endl;
    width = w;
    height = h;

    return 0;
}