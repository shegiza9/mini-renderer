# mini-renderer
图形学，软渲染器

## 项目概述
* 1、采用phong着色，支持线框，颜色，纹理模式渲染.

* 2、用tiny_obj_loader读取外部纹理，支持双线性插值，支持漫反射贴图、高光贴图和法线贴图.

* 3、实现3D裁剪，背面剔除，利用模板测试实现轮廓描边.

* 4、使用点光源实现blin-phong光照模型，并利用PCSS实现软阴影.

* 5、项目采用右手系坐标，实现了所需的数学库，使用SDL2的像素点渲染功能.
##项目阐述
render.h为主要的渲染部分，mymath.h包含了向量、矩阵、顶点、三角形等的结构体定义，main.cpp包含SDL2的使用和模型数据的导入，以及渲染实现.
用的库包括SDL2,libpng,可前往官网自行下载，项目环境为windows,visual studio2022,c++14.0.
## 效果图
* nearz top bottom left right平面裁剪
![nearz平面裁剪](https://github.com/shegiza9/images/blob/main/minirenderer/nearz%20culling.JPG)
* 漫反射+高光贴图的phong着色图
![phong](https://github.com/shegiza9/images/blob/main/minirenderer/specular%20map1.JPG)
* normal mapping(待优化，未对顶点处tangent取平均再逐像素插值）
![normal mapping](https://github.com/shegiza9/images/blob/main/minirenderer/normal%20mapping.JPG)
* shadow map深度图
![depth](https://github.com/shegiza9/images/blob/main/minirenderer/shadowmap.JPG)
* pcss软阴影
![pcss](https://github.com/shegiza9/images/blob/main/minirenderer/shadow8.JPG)
* 轮廓描边
![stencil test](https://github.com/shegiza9/images/blob/main/minirenderer/outlining2.JPG)
