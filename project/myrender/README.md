线框模式渲染：

![线框渲染](F:\GRAPHICS\output_imag\wireframe.JPG)



开启背面剔除

![背面剔除](F:\GRAPHICS\output_imag\back culling.JPG)

纹理贴图

![纹理贴图](F:\GRAPHICS\output_imag\bilinear_texture.JPG)



Gouraud着色

![gouraud着色](F:\GRAPHICS\output_imag\gouraud.JPG)



phong shading  着色更加自然，高光更明显

![phong shading](F:\GRAPHICS\output_imag\phong shading.JPG)

blin-phong shading with texture

![phong shading with texture](F:\GRAPHICS\output_imag\phong shading with texture.JPG)

由于之前简单粗暴的裁剪算法，对于一个三角形，只要有一个在视锥体外，整个三角形就被排除在渲染范围外了，所以会导致以下情况

![culling problem](F:\GRAPHICS\output_imag\culling problem.JPG)

实现对超过屏幕范围的(x,y)裁剪优化后

![裁剪优化](F:\GRAPHICS\output_imag\culling correct.JPG)







上面的图有几处白色的瑕疵，此为背面剔除时的判断方法所致，浮点数进行比较时的精度问题可能致判断失误，修正后：

![](F:\GRAPHICS\output_imag\culling 2.JPG)

近裁面问题

![](F:\GRAPHICS\output_imag\culling 3.JPG)

加入近裁面裁剪后

![](F:\GRAPHICS\output_imag\nearz culling.JPG)



![](F:\GRAPHICS\output_imag\nearz 2.JPG)

加入高光贴图后

![高光贴图1](F:\GRAPHICS\output_imag\specular map1.JPG)



未加高光贴图时的背面效果

![](F:\GRAPHICS\output_imag\back1.JPG)

高光贴图背面效果

![高光贴图2](F:\GRAPHICS\output_imag\specular map2.JPG)

法线贴图（不太自然，因为没有对顶点处的tangent和bitangent取平均，过渡不自然，但是能看出法线贴图让渲染效果更“立体”了）

![法线贴图](F:\GRAPHICS\output_imag\normal map.JPG)



加入地平面

![](F:\GRAPHICS\output_imag\with ground.JPG)

从光源处看物体得到深度图，灰度越深，深度越小，用的是点光源，但是只考虑了一个平面，没有考虑立方体深度贴图

![深度图](F:\GRAPHICS\output_imag\shadowmap.JPG)



硬阴影 走样

![走样](F:\GRAPHICS\output_imag\shadow aliasing.JPG)



考虑浮点精度后

![反走样](F:\GRAPHICS\output_imag\shadow1.JPG)



调整下光源的位置

![](F:\GRAPHICS\output_imag\shadow2.JPG)

存在的问题，锯齿，不自然

![](F:\GRAPHICS\output_imag\shadow4.JPG)

阴影应该只影响diffuse和specular，不影响ambient,修改下参数，

![](F:\GRAPHICS\output_imag\shadow6.JPG)

应用3×3的PCF（percentage-closer filtering)

![](F:\GRAPHICS\output_imag\PCF.JPG)



PCF5×5

![](F:\GRAPHICS\output_imag\PCF5.JPG)

软阴影PCSS（percentage closer soft shadows)

![](F:\GRAPHICS\output_imag\shadow7.JPG)

不太明显，换个光源位置，选光源为7×7

![](F:\GRAPHICS\output_imag\shadow8.JPG)

轮廓描边

![outlining](F:\GRAPHICS\output_imag\outlining2.JPG)

软渲染器

1、采用phong着色，支持线框，颜色，纹理模式渲染。

2、用tiny_obj_loader读取外部纹理，支持双线性插值，支持漫反射贴图、高光贴图和法线贴图

3、实现3D裁剪，背面剔除，利用模板测试实现轮廓描边

4、使用点光源实现blin-phong光照模型，并利用PCSS实现软阴影.

5、项目采用右手系坐标，使用了SDL2
