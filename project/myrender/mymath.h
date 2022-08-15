#pragma once
#include<cmath>
#include<iostream>


//设置、x不小于min，不大于max
#define PI 3.141592653
#define angle_to_radian(X) ((X)/180*PI)
#define radian_to_angle(X) ((X)/PI*180)
//计算插值：t 为 [0, 1] 之间的数值
#define interp(a, b, t)  ((1 - (t)) * (a) + (t) * (b));
typedef unsigned int UINT32;
int CMID(int x, int min, int max) { return (x < min) ? min : ((x > max) ? max : x); }
enum Render_state { WIREFRAME, TEXTURE, COLOR };


//采用右手系坐标
//vector类
//4维向量

class Vec4f {
public:
    Vec4f() = default;
    Vec4f(float _x, float _y, float _z, float _w) :x(_x), y(_y), z(_z), w(_w) {};

    friend Vec4f operator+(const Vec4f& a, const Vec4f& b);
    friend Vec4f operator-(const Vec4f& a, const Vec4f& b);

    Vec4f& operator+=(const Vec4f& v2);
    Vec4f& operator-=(const Vec4f& v2);


    float length()const {
        return (float)sqrt(x * x + y * y + z * z);
    }

    float square_length()const {
        return x * x + y * y + z * z;
    }

    float dotProduct(const Vec4f& a) const {
        return x * a.x + y * a.y + z * a.z;
    }

    Vec4f crossProduct(const Vec4f& a)const {
        return Vec4f(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x, 1.0f);
    }
    //归一化
    void normalized() {
        float len = this->length();
        if (len != 0.0f) {
            float invlen = 1 / len;
            x *= invlen;
            y *= invlen;
            z *= invlen;
        }
    }
    Vec4f normalize() {
        float len = this->length();
        Vec4f ret = { 0,0,0,1 };
        if (len != 0.0f) {
            float invlen = 1 / len;
            ret = { x * invlen, y * invlen, z * invlen, 1.0f };
            return ret;
        }
        return ret;
    }

    //插值(1-t)*this+t*v
    Vec4f vec_interpolate(const Vec4f& v, float t)const {
        Vec4f ret;
        ret.x = interp(x, v.x, t);
        ret.y = interp(y, v.y, t);
        ret.z = interp(z, v.z, t);
        ret.w = 1.0f;
        return ret;
    }
    Vec4f operator*(float k)const {
        Vec4f ret(x * k, y * k, z * k, 1.0f);
        return ret;
    }
    Vec4f operator*(const Vec4f& a) {
        Vec4f ret(x * a.x, y * a.y, z * a.z, 1.0f);
        return ret;
    }
    float x, y, z, w;
    ~Vec4f() = default;
};
inline Vec4f& Vec4f::operator+=(const Vec4f& v2) {
    x += v2.x;
    y += v2.y;
    z += v2.z;
    w = 1.0f;
    return *this;
}
inline Vec4f& Vec4f::operator-=(const Vec4f& v2) {
    x -= v2.x;
    y -= v2.y;
    z -= v2.z;
    w = 1.0f;
    return *this;
}
inline Vec4f operator+(const Vec4f& a, const Vec4f& b) {
    return Vec4f(a.x + b.x, a.y + b.y, a.z + b.z, 1.0f);
}
inline Vec4f operator-(const Vec4f& a, const Vec4f& b) {
    return Vec4f(a.x - b.x, a.y - b.y, a.z - b.z, 1.0f);
}



//四维矩阵
struct Matrix {
	Matrix() = default;
	Matrix operator+(const Matrix&)const;
	Matrix operator-(const Matrix&)const;
	Matrix operator*(const Matrix&)const;
	Matrix operator*(float k)const;
	Matrix transpose()const;
	Matrix inverse3f()const;
	Matrix inverse4f()const;
	Vec4f operator*(const Vec4f&)const;
	float m[4][4];
};

Matrix Matrix::operator+(const Matrix& a)const {
	Matrix b;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			b.m[i][j] = this->m[i][j] + a.m[i][j];
		}
	}
	return b;
}

Matrix Matrix::operator-(const Matrix& a)const {
	Matrix b;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			b.m[i][j] = this->m[i][j] - a.m[i][j];
		}
	}
	return b;
}

Matrix Matrix::operator*(const Matrix& a)const {
	Matrix b;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			b.m[i][j] = m[i][0] * a.m[0][j] +
				m[i][1] * a.m[1][j] +
				m[i][2] * a.m[2][j] +
				m[i][3] * a.m[3][j];
		}
	}
	return b;
}

Matrix Matrix::operator*(float k)const {
	Matrix a;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			a.m[i][j] = m[i][j] * k;
		}
	}
	return a;
}


Vec4f Matrix::operator*(const Vec4f& x)const {
	Vec4f y;
	float X = x.x, Y = x.y, Z = x.z, W = x.w;
	y.x = X * m[0][0] + Y * m[0][1] + Z * m[0][2] + W * m[0][3];
	y.y = X * m[1][0] + Y * m[1][1] + Z * m[1][2] + W * m[1][3];
	y.z = X * m[2][0] + Y * m[2][1] + Z * m[2][2] + W * m[2][3];
	y.w = X * m[3][0] + Y * m[3][1] + Z * m[3][2] + W * m[3][3];
	return y;
}

Matrix Matrix::transpose()const {
	Matrix ret;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			ret.m[i][j] = this->m[j][i];
		}
	}
	return ret;
}

Matrix Matrix::inverse3f()const {
	Matrix ret;
	float t[3][6];
	int i, j, k;
	float f;
	//高斯消元法
	//对左上3乘3矩阵进行求逆
	for (i = 0; i < 3; i++)
		for (j = 0; j < 6; j++) {
			if (j < 3)
				t[i][j] = m[i][j];
			else if (j == i + 3)
				t[i][j] = 1;
			else
				t[i][j] = 0;
		}
	//把左边变成单位矩阵
	for (i = 0; i < 3; i++) {
		f = t[i][i];
		for (j = 0; j < 6; j++)
			t[i][j] /= f;
		for (j = 0; j < 3; j++) {
			if (j != i) {
				f = t[j][i];
				for (k = 0; k < 6; k++)
					t[j][k] = t[j][k] - t[i][k] * f;
			}
		}
	}

	for (i = 0; i < 3; i++)
		for (j = 3; j < 6; j++)
			ret.m[i][j - 3] = t[i][j];

	/*ret.m[3][0] = -m[3][0];
	ret.m[3][1] = -m[3][1];
	ret.m[3][2] = -m[3][2];*/

	ret.m[3][0] = ret.m[3][1] = ret.m[3][2] = 0;
	ret.m[0][3] = ret.m[1][3] = ret.m[2][3] = 0;
	ret.m[3][3] = 1.0f;

	return ret;
}

Matrix Matrix::inverse4f()const {
	Matrix ret;
	float t[4][8];
	int i, j, k;
	float f;

	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 8; ++j) {
			if (j < 4) {
				t[i][j] = m[i][j];
			}
			else if (j == i + 4) {
				t[i][j] = 1.0f;
			}
			else {
				t[i][j] = 0;
			}
		}
	}

	for (i = 0; i < 4; ++i) {
		f = t[i][i];
		for (j = 0; j < 8; ++j) {
			t[i][j] /= f;
		}
		for (j = 0; j < 4; ++j) {
			if (j != i) {
				f = t[j][i];
				for (k = 0; k < 8; ++k) {
					t[j][k] = t[j][k] - t[i][k] * f;
				}
			}
		}
	}

	for (i = 0; i < 4; ++i) {
		for (j = 4; j < 8; ++j) {
			ret.m[i][j - 4] = t[i][j];
		}
	}
	return ret;
}
//设置零矩阵
Matrix matrix_set_zero() {
	Matrix m;
	m.m[0][0] = m.m[0][1] = m.m[0][2] = m.m[0][3] = 0.0f;
	m.m[1][0] = m.m[1][1] = m.m[1][2] = m.m[1][3] = 0.0f;
	m.m[2][0] = m.m[2][1] = m.m[2][2] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = m.m[3][3] = 0.0f;

	return m;
}
//设置单位矩阵
Matrix matrix_set_identity() {
	Matrix m;
	m.m[0][0] = m.m[1][1] = m.m[2][2] = m.m[3][3] = 1.0f;
	m.m[0][1] = m.m[0][2] = m.m[0][3] = 0.0f;
	m.m[1][0] = m.m[1][2] = m.m[1][3] = 0.0f;
	m.m[2][0] = m.m[2][1] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;

	return m;
}
//平移变换
Matrix matrix_set_translate(float a, float b, float c) {
	Matrix m = matrix_set_identity();
	m.m[0][3] = a;
	m.m[1][3] = b;
	m.m[2][3] = c;

	return m;
}
//缩放变化
Matrix matrix_set_scale(float a, float b, float c) {
	Matrix m = matrix_set_identity();
	m.m[0][0] = a;
	m.m[1][1] = b;
	m.m[2][2] = c;

	return m;
}
//旋转矩阵，用四元数
Matrix matrix_set_rotate(float x, float y, float z, float theta) {
	Matrix m;
	float qsin = (float)sin(theta * 0.5f);
	float qcos = (float)cos(theta * 0.5f);
	Vec4f vec(x, y, z, 1.0f);
	vec.normalized();
	x = vec.x * qsin;
	y = vec.y * qsin;
	z = vec.z * qsin;
	m.m[0][0] = 1 - 2 * y * y - 2 * z * z;
	m.m[0][1] = 2 * x * y - 2 * z * qcos;
	m.m[0][2] = 2 * x * z + 2 * y * qcos;
	m.m[1][0] = 2 * x * y + 2 * z * qcos;
	m.m[1][1] = 1 - 2 * x * x - 2 * z * z;
	m.m[1][2] = 2 * y * z - 2 * x * qcos;
	m.m[2][0] = 2 * x * z - 2 * y * qcos;
	m.m[2][1] = 2 * y * z + 2 * x * qcos;
	m.m[2][2] = 1 - 2 * x * x - 2 * y * y;
	m.m[0][3] = m.m[1][3] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
	m.m[3][3] = 1.0f;

	return m;
}
//设置摄像机矩阵,右手系
Matrix matrix_set_lookat(const Vec4f& eye, const Vec4f& at, const Vec4f& up) {
	Matrix m;
	Vec4f xaxis, yaxis, zaxis;
	zaxis = at;
	zaxis.normalized();
	xaxis = up.crossProduct(zaxis);
	xaxis.normalized();
	yaxis = zaxis.crossProduct(xaxis);
	m.m[0][0] = xaxis.x;
	m.m[0][1] = xaxis.y;
	m.m[0][2] = xaxis.z;
	m.m[0][3] = -xaxis.dotProduct(eye);

	m.m[1][0] = yaxis.x;
	m.m[1][1] = yaxis.y;
	m.m[1][2] = yaxis.z;
	m.m[1][3] = -yaxis.dotProduct(eye);

	m.m[2][0] = zaxis.x;
	m.m[2][1] = zaxis.y;
	m.m[2][2] = zaxis.z;
	m.m[2][3] = -zaxis.dotProduct(eye);

	m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
	m.m[3][3] = 1.0f;

	return m;
}
//设置透视投影矩阵,nz和fz均为正数，需要转换为负数.
Matrix matris_set_perspective(float fov, float aspect, float nz, float fz) {
	Matrix m;
	m = matrix_set_identity();
	float inv_tanfov = 1 / tan(fov / 2);
	float inv_asp = 1 / aspect;

	m.m[0][0] = inv_asp * inv_tanfov;
	m.m[1][1] = inv_tanfov;
	m.m[2][2] = -(nz + fz) / (fz - nz);
	m.m[2][3] = -2 * fz * nz / (fz - nz);
	m.m[3][2] = -1;
	m.m[3][3] = 0.0f;
	return m;
}


//投影变换矩阵类

class TransformMatrix {
public:
	TransformMatrix() {
		world = matrix_set_identity();
		view = matrix_set_identity();
		perspective = matris_set_perspective(PI * 0.5f, 1.333f, 1.0f, 500.0f);
		viewdirection = { 1.0f,0.0f,0.0f,0.0f };
		width = 800;
		height = 600;
		eye = { 0,0,0,1.0f };
		nz = 1.0f;
		fz = 500.0f;
		tranform_update();
	}
	void set_world_rotate(float x, float y, float z, float theta) {
		world = matrix_set_rotate(x, y, z, theta);
	}
	void tranform_update() {
		m_matrix = world;
		vm_matrix = view * world;
		pv_matrix = perspective * view;
		transform = pv_matrix * world;
		view_inv = view.inverse4f();
		vm_inv_tran = vm_matrix.inverse3f().transpose();
	}
	TransformMatrix(int w, int h) {
		world = matrix_set_identity();
		view = matrix_set_identity();
		float aspect = (float)w / (float)h;
		perspective = matris_set_perspective(PI * 0.5f, 1.333f, 1.0f, 500.0f);
		width = w;
		height = h;
		tranform_update();
	}
	TransformMatrix(const TransformMatrix& trans) {
		world = trans.world;
		view = trans.view;
		perspective = trans.perspective;
		transform = trans.transform;
		width = trans.width;
		height = trans.height;
	}
	TransformMatrix operator=(const TransformMatrix& trans) {
		world = trans.world;
		view = trans.view;
		perspective = trans.perspective;
		transform = trans.transform;
		width = trans.width;
		height = trans.height;
		return *this;
	}
	void trans_apply(Vec4f& y, const Vec4f& x) {
		y = transform * x;
	}
	//得到屏幕坐标
	void transform_homogenize(Vec4f& y, const Vec4f& x) {
		float inv_w = 1 / x.w;
		y.x = (x.x * inv_w + 1.0f) * 0.5 * width;
		y.y = (1.0f - x.y * inv_w) * 0.5 * height;
		y.z = x.z * inv_w;
		y.w = 1.0f;
	}

	void set_viewdir(const Vec4f& eye, const Vec4f& at) {
		this->eye = eye;

		viewdirection = at;
	}

	int transform_check_cvv(const Vec4f& v) {
		float w = v.w;
		int check = 0;

		if (v.z < -w) check |= 1;
		if (v.z > w) check |= 2;
		if (v.x < -w) check |= 4;
		if (v.x > w)check |= 8;
		if (v.y < -w) check |= 16;
		if (v.y > w) check |= 32;
		return check;
	}
	Matrix world;
	Matrix view;
	Matrix perspective;
	Vec4f eye;

	Vec4f viewdirection;
	float nz, fz;
	float width, height;
	Matrix m_matrix;
	Matrix vm_matrix;
	Matrix pv_matrix;
	Matrix transform;//transform=perspective*view*world;
	Matrix vm_inv_tran;
	Matrix view_inv;
};


//RGB颜色
class Color {
public:
	Color() = default;
	Color(float rr, float gg, float bb) :r(rr), g(gg), b(bb) {};
	Color operator*(const Vec4f& k) {
		Color ret(k.x * r, k.y * g, k.z * b);
		return ret;
	}
	Color operator*(float s) {
		Color ret(r * s, g * s, b * s);
		return ret;
	}

	friend Color operator+(const Color& c1, const Color& c2) {
		Color ret;
		ret.r = c1.r + c2.r;
		ret.g = c1.g + c2.g;
		ret.b = c1.b + c2.b;
		return ret;
	}
	float r, g, b;
};
typedef Vec4f Pos;
//纹理坐标
struct Texturecoord {
	float u, v;
};
//顶点：位置坐标，纹理坐标，颜色，透视矫正
struct Vertex_t {
	Pos pos;
	Texturecoord tc;
	Color color;
	Vec4f normal;
	float rhw;
};
//透视矫正
void vertex_init_rhw(Vertex_t& v) {
	float rhw = 1.0f / v.pos.w;
	v.rhw = rhw;
	/*v.tc.u *= rhw;
	v.tc.v *= rhw;
	v.color.r *= rhw;
	v.color.g *= rhw;
	v.color.b *= rhw;*/
}
//透视矫正插值，input的ver1,ver2都是经过归一化的（transform_homogenize），
//且ver(1/2).pos.w保留着原来的w,然后经过了vertex_init_rhw函数处理
void vertex_interp(Vertex_t& y, const Vertex_t& ver1, const Vertex_t& ver2, float t) {
	y.pos = ver1.pos.vec_interpolate(ver2.pos, t);
	y.rhw = interp(ver1.rhw, ver2.rhw, t);
	/*y.tc.u = interp(ver1.tc.u, ver2.tc.u, t);
	y.tc.v = interp(ver1.tc.v, ver2.tc.v, t);
	y.color.r = interp(ver1.color.r, ver2.color.r, t);
	y.color.g = interp(ver1.color.g, ver2.color.g, t);
	y.color.b = interp(ver1.color.b, ver2.color.b, t);*/
}
//用于扫描线绘制，分割扫描线
void vertex_division(Vertex_t& y, const Vertex_t& ver1, const Vertex_t& ver2, float w) {
	float inv_w = 1.0f / w;
	y.pos = (ver2.pos - ver1.pos) * inv_w;
	/*y.tc.u = (ver2.tc.u -ver1.tc.u) * inv_w;
	y.tc.v = (ver2.tc.v- ver1.tc.v) * inv_w;
	y.color.r = (ver2.color.r - ver1.color.r) * inv_w;
	y.color.g = (ver2.color.g - ver1.color.g) * inv_w;
	y.color.b = (ver2.color.b - ver1.color.b) * inv_w;*/
	y.rhw = (ver2.rhw - ver1.rhw) * inv_w;
}
//逐像素相加
void vertex_add(Vertex_t& y, const Vertex_t& x) {
	y.pos += x.pos;
	y.rhw += x.rhw;
	/*y.tc.u += x.tc.u;
	y.tc.v += x.tc.v;
	y.color.r += x.color.r;
	y.color.g += x.color.g;
	y.color.b += x.color.b;*/
}

struct Edge { Vertex_t v, v1, v2; };
struct Trapezoid { float top, bottom; Edge left, right; };
//v为起始点（为float类型），x,y为其像素坐标，w为扫描线像素长度，step为每走一步所要加的vertex属性
struct Scanline { Vertex_t v, step; int x, y, w; };
struct Screen_points { Pos p1, p2, p3; };

class Triangle {
public:
	Triangle() = default;
	void set_vertex_pos(int i, float x, float y, float z) {
		this->v[i].pos = { x,y,z,1.0f };
	}
	void set_vertex_texcord(int i, float u, float v) {
		this->v[i].tc = { u,v };
	}
	void set_vertex_color(int i, float r, float g, float b) {
		this->v[i].color = { r,g,b };
	}
	void set_vertex_rhw(int i, float rhw) {
		this->v[i].rhw = rhw;
	}
	void set_vertex_normal(int i, float nx, float ny, float nz) {
		this->v[i].normal = { nx,ny,nz,1.0f };
	}

	~Triangle() = default;
	Vertex_t v[3];

};

/****扫描法画三角形***/
//分割
int trapezoid_init_triangle(Trapezoid* trap, const Vertex_t* p1, const Vertex_t* p2, const Vertex_t* p3) {
	const Vertex_t* p;
	//p1,p2,p3,y值从小到大
	if (p1->pos.y > p2->pos.y) { p = p1; p1 = p2; p2 = p; }
	if (p1->pos.y > p3->pos.y) { p = p1; p1 = p3; p3 = p; }
	if (p2->pos.y > p3->pos.y) { p = p2; p2 = p3; p3 = p; }
	if (p1->pos.y == p2->pos.y && p1->pos.y == p3->pos.y) return 0;
	if (p1->pos.x == p2->pos.x && p1->pos.x == p3->pos.x) return 0;

	//triangle down
	if (p1->pos.y == p2->pos.y) {
		if (p1->pos.x > p2->pos.x) p = p1, p1 = p2, p2 = p;
		trap[0].top = p1->pos.y;
		trap[0].bottom = p3->pos.y;
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p3;
		trap[0].right.v1 = *p2;
		trap[0].right.v2 = *p3;
		return (trap[0].top < trap[0].bottom) ? 1 : 0;
	}
	//triangle up
	if (p2->pos.y == p3->pos.y) {
		if (p2->pos.x > p3->pos.x) p = p2, p2 = p3, p3 = p;
		trap[0].top = p1->pos.y;
		trap[0].bottom = p3->pos.y;
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p2;
		trap[0].right.v1 = *p1;
		trap[0].right.v2 = *p3;
		return (trap[0].top < trap[0].bottom) ? 1 : 0;
	}
	trap[0].top = p1->pos.y;
	trap[0].bottom = p2->pos.y;
	trap[1].top = p2->pos.y;
	trap[1].bottom = p3->pos.y;

	float k = (p3->pos.y - p1->pos.y) / (p2->pos.y - p1->pos.y);
	float x = p1->pos.x + (p2->pos.x - p1->pos.x) * k;

	//p2在最左边
	if (x <= p3->pos.x) {
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p2;
		trap[0].right.v1 = *p1;
		trap[0].right.v2 = *p3;

		trap[1].left.v1 = *p2;
		trap[1].left.v2 = *p3;
		trap[1].right.v1 = *p1;
		trap[1].right.v2 = *p3;
	}
	//p2在最右边
	else {
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p3;
		trap[0].right.v1 = *p1;
		trap[0].right.v2 = *p2;

		trap[1].left.v1 = *p1;
		trap[1].left.v2 = *p3;
		trap[1].right.v1 = *p2;
		trap[1].right.v2 = *p3;
	}
	return 2;
}
void trapezoid_edge_interp(Trapezoid& trap, float y) {
	float s1 = trap.left.v2.pos.y - trap.left.v1.pos.y;
	float s2 = trap.right.v2.pos.y - trap.right.v1.pos.y;
	float t1 = (y - trap.left.v1.pos.y) / s1;
	float t2 = (y - trap.right.v1.pos.y) / s2;
	vertex_interp(trap.left.v, trap.left.v1, trap.left.v2, t1);
	vertex_interp(trap.right.v, trap.right.v1, trap.right.v2, t2);
}
void trapezoid_init_scan_line(const Trapezoid& trap, Scanline& scanline, int y) {
	float width = trap.right.v.pos.x - trap.left.v.pos.x;
	scanline.x = (int)(trap.left.v.pos.x + 0.5f);
	scanline.w = (int)(trap.right.v.pos.x + 0.5f) - scanline.x;
	scanline.y = y;
	scanline.v = trap.left.v;
	if (trap.left.v.pos.x >= trap.right.v.pos.x) scanline.w = 0;
	vertex_division(scanline.step, trap.left.v, trap.right.v, width);
}
/****扫描法画三角形***/