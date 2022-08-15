# graphics-notes
图形学，软渲染器

## 项目概述
1、采用phong着色，支持线框，颜色，纹理模式渲染.

2、用tiny_obj_loader读取外部纹理，支持双线性插值，支持漫反射贴图、高光贴图和法线贴图.

3、实现3D裁剪，背面剔除，利用模板测试实现轮廓描边.

4、使用点光源实现blin-phong光照模型，并利用PCSS实现软阴影.

5、项目采用右手系坐标，实现了所需的数学库，使用SDL2的像素点渲染功能.
![背面剔除](https://github.com/shegiza9/images/raw/main/minirenderer/back%20culling.JPG)
![nearz平面裁剪](https://github.com/shegiza9/images/blob/main/minirenderer/nearz%20culling.JPG)
