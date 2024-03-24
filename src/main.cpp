#define _USE_MATH_DEFINES


#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>
using Eigen::VectorXd;


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STBI_MSC_SECURE_CRT
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"




class material {
public:
	float spec[3] = { 1, 1, 1 };
	float diff[3] = { 1, 0, 0 };
	float shiny = 200;

	material() {}

	material(float* diff, float* spec) {
		this->diff[0] = diff[0];
		this->diff[1] = diff[1];
		this->diff[2] = diff[2];
		this->spec[0] = spec[0];
		this->spec[1] = spec[1];
		this->spec[2] = spec[2];
	}
};

class sphere {
public:
	VectorXd o;
	int r;
	material mat = material();

	sphere() {}
	sphere(VectorXd o, float r) {
		this->o = o;
		this->r = r;
		this->mat = material();
	}
	sphere(float x, float y, float z, float w, float r) {
		VectorXd o(x, y, z, w);
		this->o = o;
		this->r = r;
		this->mat = material();
	}
	sphere(float x, float y, float z, float w, float r, material mat) {
		VectorXd o(x, y, z, w);
		this->o = o;
		this->r = r;
		this->mat = mat;
	}
};

class hyperplane {
public:
	VectorXd p1;
	VectorXd p2;
	VectorXd p3;
	VectorXd p4;
	VectorXd n;
	material mat = material();

	hyperplane() {}
	hyperplane(VectorXd p1, VectorXd p2, VectorXd p3, VectorXd p4) {
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
		this->p4 = p4;
		VectorXd u = p2 - p1;
		VectorXd v = p3 - p1;
		VectorXd w = p4 - p1;
		this->n = hypernormal(u, v, w);
	}
	hyperplane(VectorXd p1, VectorXd p2, VectorXd p3, VectorXd p4, material mat) {
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
		this->p4 = p4;
		VectorXd u = p2 - p1;
		VectorXd v = p3 - p1;
		VectorXd w = p4 - p1;
		this->n = hypernormal(u, v, w);
		this->mat = mat;
	}
};

VectorXd hypernormal(VectorXd u, VectorXd v, VectorXd w) {

};

class light {
public:
	VectorXd o;
	material mat = material();

	light() {}
	light(VectorXd o, material mat) {
		this->o = o;
		this->mat = mat;
	}
	light(VectorXd o) {
		this->o = o;
	}
	light(VectorXd o, float* diff) {
		this->o = o;
		float spec[3] = { 1, 1, 1 };
		this->mat = material(diff, spec);
	}
};

std::vector<sphere> spheres;
std::vector<light> lights;
std::vector<hyperplane> triangles;
float amb = 2;
float rr = 0.5;
float bg[3] = { 1,0.2,1 };
float maxCol = 0;

class ray3 {
public:
	VectorXd v;
	VectorXd o;

	ray3(VectorXd v, VectorXd o) {
		this->v = v;
		this->o = o;
	}
	ray3() {}
};


float sphereIntersect(ray3 ray, sphere s, bool debug = false) {
	VectorXd v;
	VectorXd o;
	float a = v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
	float b = 2 * (v.x() * (o.x() - s.o.x()) + v.y() * (o.y() - s.o.y()) + v.z() * (o.z() - s.o.z()));
	float c = (o.x() - s.o.x()) * (o.x() - s.o.x()) + (o.y() - s.o.y()) * (o.y() - s.o.y()) + (o.z() - s.o.z()) * (o.z() - s.o.z()) - (s.r * s.r);
	float d = b * b - 4 * a * c;
	//printf("\nv, o: (%f,%f,%f), (%f, %f, %f)\n", v.x, v.y, v.z, o.x, o.y, o.z);
	//printf("a, b, c, d: %f, %f, %f, %f\n", a, b, c, d);
	if (d < 0) {
		if (debug) printf("Not on ray's path\n");
		return NULL;
	}
	else if (d == 0) {
		if (debug) printf("One intersection\n");
		if (-b < 0) return NULL;
		else return -b;
	}
	else {
		float t1 = (-b + sqrt(d)) / (2 * a);
		float t2 = (-b - sqrt(d)) / (2 * a);
		if (t1 < 0 && t2 < 0) {
			if (debug) printf("Sphere behind\n");
			return NULL;
		}
		else if ((t1 < 0 && t2 > 0) || (t2 < 0 && t1 > 0)) {
			if (debug) printf("Inside sphere\n");
			return -1;
		}
		else {
			//printf("\nray: (%f, %f, %f) -> (%f, %f, %f)\n", ray.o.x, ray.o.y, ray.o.z, ray.v.x, ray.v.y, ray.v.z);
			//printf("t1: %f    t2: %f\n", t1, t2);
			if (t1 < t2) {
				//printf("to_o: %f t1: %f\n", sqrt(pow(s.o.x - ray.o.x, 2) + pow(s.o.y - ray.o.y, 2) + pow(s.o.z - ray.o.z, 2)), t1);
				return t1;
			}
			else {
				//printf("to_o: %f t2: %f\n", sqrt(pow(s.o.x - ray.o.x, 2) + pow(s.o.y - ray.o.y, 2) + pow(s.o.z - ray.o.z, 2)), t2);
				return t1;
			}
		}
	}
}



//! Paused updating here


float hyperplaneIntersect(ray3 ray, hyperplane tri) {
	float temp = glm::dot(tri.n, ray.v);
	if (temp == 0) return NULL;
	float d = glm::dot(tri.n, tri.p1);
	float t = (d - glm::dot(tri.n, ray.o)) / temp;
	glm::vec3 pt = ray.o + (ray.v * t);

	float n1 = glm::dot(glm::cross(tri.p2 - tri.p1, pt - tri.p1), tri.n);
	float n2 = glm::dot(glm::cross(tri.p3 - tri.p2, pt - tri.p2), tri.n);
	float n3 = glm::dot(glm::cross(tri.p1 - tri.p3, pt - tri.p3), tri.n);
	if (n1 < 0 || n2 < 0 || n3 < 0) {
		return NULL;
	}
	else {
		return t;
	}
}

ray3 reflect(ray3 r, glm::vec3 n, glm::vec3 pt) {
	glm::vec3 dir = glm::reflect(r.v, n);
	return ray3(dir, pt);
}

float* traceRay(ray3 ray, int depth) {
	material mat = material();
	float dist = INFINITY;
	bool hit = false;
	glm::vec3 normal;
	glm::vec3 pt;
	sphere sp;
	triangle tr;

	for (sphere s : spheres) {
		float t = sphereIntersect(ray, s);
		if (t > 0 && t < dist) {
			hit = true;
			sp = s;
			pt = glm::vec3(ray.o.x + ray.v.x * t, ray.o.y + ray.v.y * t, ray.o.z + ray.v.z * t);
			normal = glm::normalize(glm::vec3(pt.x - s.o.x, pt.y - s.o.y, pt.z - s.o.z));
			pt.x += normal.x * 0.0001;
			pt.y += normal.y * 0.0001;
			pt.z += normal.z * 0.0001;
			dist = t;
			mat = s.mat;
		}
	}
	for (triangle tri : triangles) {
		float t = triangleIntersect(ray, tri);
		if (t > 0 && t < dist) {
			sp = sphere(1, 1, 1, -1);
			tr = tri;
			hit = true;
			pt = glm::vec3(ray.o.x + ray.v.x * t, ray.o.y + ray.v.y * t, ray.o.z + ray.v.z * t);
			normal = -tri.n;
			pt.x += normal.x * 0.0001;
			pt.y += normal.y * 0.0001;
			pt.z += normal.z * 0.0001;
			dist = t;
			mat = tri.mat;
		}
	}


	if (dist != INFINITY) {
		ray3 reflectedRay;
		reflectedRay = reflect(ray, normal, pt);
		float color[3] = { mat.diff[0] * amb, mat.diff[1] * amb, mat.diff[2] * amb };

		for (light l : lights) {
			float distToL = glm::distance(pt, l.o);
			ray3 r(glm::normalize(glm::vec3(l.o.x - pt.x, l.o.y - pt.y, l.o.z - pt.z)), pt);
			bool lit = true;
			for (sphere s : spheres) {
				if (s.o == sp.o) continue;
				float tempDist = sphereIntersect(r, s);
				//printf("distToL: %f   tempDist: %f\n", distToL, tempDist);
				if (tempDist && tempDist < distToL) {
					lit = false;
					//printf("Obstructed\n");
					break;
				}
				//else printf("Not obstructing\n");
			}
			if (lit) {
				//printf("lit: (%f, %f, %f)\n", color[0], color[1], color[2]);
				float cont = pow(fmax(glm::dot(glm::normalize(r.v + -ray.v), normal), 0), mat.shiny);
				//if (cont > 0.1) printf("Cont: %f\n", cont);
				color[0] += mat.spec[0] * l.mat.diff[0] * cont * 30 + mat.diff[0] * amb;
				color[1] += mat.spec[1] * l.mat.diff[1] * cont * 30 + mat.diff[1] * amb;
				color[2] += mat.spec[2] * l.mat.diff[2] * cont * 30 + mat.diff[2] * amb;
				//printf("after: (%f, %f, %f)\n\n", color[0], color[1], color[2]);
			}
		}


		if (depth <= 0 || (mat.diff[0] == 1 && mat.diff[1] == 1 && mat.diff[2] == 1)) {
			if (color[0] > maxCol) maxCol = color[0];
			if (color[1] > maxCol) maxCol = color[1];
			if (color[2] > maxCol) maxCol = color[2];
			return color;
		}


		float* newCol = traceRay(reflectedRay, depth - 1);

		color[0] += newCol[0] * rr;
		color[1] += newCol[1] * rr;
		color[2] += newCol[2] * rr;

		if (color[0] > maxCol) maxCol = color[0];
		if (color[1] > maxCol) maxCol = color[1];
		if (color[2] > maxCol) maxCol = color[2];



		return color;
	}
	else {
		return bg;
	}
}


const int outWidth = 800;
const int outHeight = 800;
float fovx = M_PI / 2;
float fovy = fovx * outHeight / outWidth;
int maxDepth = 4;

int main(int argc, char* argv[]) {

	std::ifstream infile;
	infile.open("scene.txt");

	if (infile.is_open()) {
		printf("Reading file. . .\n");
		std::string str;
		while (std::getline(infile, str)) {
			if (str.size() > 0 && str.at(0) == '-') {
				std::string temp;
				if (str == "-BACKGROUND") {
					printf("Background found\n");
					// Background Color
					std::getline(infile, temp);
					bg[0] = stof(temp);
					std::getline(infile, temp);
					bg[1] = stof(temp);
					std::getline(infile, temp);
					bg[2] = stof(temp);
				}
				else if (str == "-MAXDEPTH") {
					printf("Max depth found\n");
					std::getline(infile, temp);
					printf("Depth read: %d\n", stoi(temp));
					maxDepth = stoi(temp);
				}
				else if (str == "-LIGHT") {
					printf("Light found\n");
					light nov = light();
					// Origin
					std::getline(infile, temp);
					nov.o.x = stof(temp);
					std::getline(infile, temp);
					nov.o.y = stof(temp);
					std::getline(infile, temp);
					nov.o.z = stof(temp);
					std::getline(infile, temp);

					// Diffuse
					std::getline(infile, temp);
					nov.mat.diff[0] = stof(temp);
					std::getline(infile, temp);
					nov.mat.diff[1] = stof(temp);
					std::getline(infile, temp);
					nov.mat.diff[2] = stof(temp);
					std::getline(infile, temp);

					// Specular
					std::getline(infile, temp);
					nov.mat.spec[0] = stof(temp);
					std::getline(infile, temp);
					nov.mat.spec[1] = stof(temp);
					std::getline(infile, temp);
					nov.mat.spec[2] = stof(temp);
					lights.push_back(nov);
				}
				else if (str == "-SPHERE") {
					printf("Sphere found\n");
					sphere nov = sphere();

					// Origin and radius
					std::getline(infile, temp);
					nov.o.x = stof(temp);
					std::getline(infile, temp);
					nov.o.y = stof(temp);
					std::getline(infile, temp);
					nov.o.z = stof(temp);
					std::getline(infile, temp);
					nov.r = stof(temp);
					std::getline(infile, temp);

					// Diffuse
					std::getline(infile, temp);
					nov.mat.diff[0] = stof(temp);
					std::getline(infile, temp);
					nov.mat.diff[1] = stof(temp);
					std::getline(infile, temp);
					nov.mat.diff[2] = stof(temp);
					std::getline(infile, temp);

					// Specular
					std::getline(infile, temp);
					nov.mat.spec[0] = stof(temp);
					std::getline(infile, temp);
					nov.mat.spec[1] = stof(temp);
					std::getline(infile, temp);
					nov.mat.spec[2] = stof(temp);
					std::getline(infile, temp);

					// Shininess
					std::getline(infile, temp);
					nov.mat.shiny = stof(temp);
					spheres.push_back(nov);
				}
				else if (str == "-QUAD") {
					printf("Quad found\n");
					triangle nov = triangle();
					triangle nov2 = triangle();

					// First point
					std::getline(infile, temp);
					nov.p1.x = stof(temp);
					std::getline(infile, temp);
					nov.p1.y = stof(temp);
					std::getline(infile, temp);
					nov.p1.z = stof(temp);
					std::getline(infile, temp);

					// Second point
					std::getline(infile, temp);
					nov.p2.x = stof(temp);
					std::getline(infile, temp);
					nov.p2.y = stof(temp);
					std::getline(infile, temp);
					nov.p2.z = stof(temp);
					std::getline(infile, temp);

					// Third point
					std::getline(infile, temp);
					nov.p3.x = stof(temp);
					std::getline(infile, temp);
					nov.p3.y = stof(temp);
					std::getline(infile, temp);
					nov.p3.z = stof(temp);
					std::getline(infile, temp);

					// Copy redundant points
					nov2.p1 = nov.p2;
					nov2.p3 = nov.p3;

					// Fourth point
					std::getline(infile, temp);
					nov2.p2.x = stof(temp);
					std::getline(infile, temp);
					nov2.p2.y = stof(temp);
					std::getline(infile, temp);
					nov2.p2.z = stof(temp);
					std::getline(infile, temp);

					// Diffuse
					std::getline(infile, temp);
					nov.mat.diff[0] = stof(temp);
					std::getline(infile, temp);
					nov.mat.diff[1] = stof(temp);
					std::getline(infile, temp);
					nov.mat.diff[2] = stof(temp);
					std::getline(infile, temp);

					// Specular
					std::getline(infile, temp);
					nov.mat.spec[0] = stof(temp);
					std::getline(infile, temp);
					nov.mat.spec[1] = stof(temp);
					std::getline(infile, temp);
					nov.mat.spec[2] = stof(temp);
					std::getline(infile, temp);

					// Shininess
					std::getline(infile, temp);
					nov.mat.shiny = stof(temp);
					nov2.mat = nov.mat;

					glm::vec3 u = nov.p2 - nov.p1;
					glm::vec3 v = nov.p3 - nov.p1;
					nov.n = glm::normalize(glm::cross(u, v));
					u = nov2.p2 - nov2.p1;
					v = nov2.p3 - nov2.p1;
					nov2.n = glm::normalize(glm::cross(u, v));

					triangles.push_back(nov);
					triangles.push_back(nov2);
				}
			}
		}
		printf("Done reading file.\n");
	}

	fovy = fovx * outHeight / outWidth;
	const int channels = 3;
	const float a = (float)outWidth / outHeight;
	float d = 1 / tan(fovy / 2);


	auto out = new unsigned char[outWidth][outHeight][3];
	glm::vec3 eye(0, 0, 0);
	glm::vec3 l(0, 0, 1);
	glm::vec3 u(0, 1, 0);
	glm::vec3 v(1, 0, 0);
	glm::vec3 ll = eye + (d * l) - (a * v) - u;
	//printf("ll: (%f,%f,%f, %f)", ll.x, ll.y, ll.z, ll.w);


	for (int i = 0; i < outWidth; i++) {
		for (int j = 0; j < outHeight; j++) {
			glm::vec3 r = 2 * a * v * ((float)i / outWidth);
			glm::vec3 d = (float)2 * u * ((float)j / outHeight);
			glm::vec3 p = ll + r + d;
			//printf("\n(%f, %f, %f) + (%f, %f, %f) = (%f, %f, %f)", r.x, r.y, r.z, d.x, d.y, d.z, p.x, p.y, p.z);
			p = glm::normalize(p);
			ray3 ray(p, eye);

			//printf("dx, dy, dz: %f, %f, %f\n", ray[0], ray[1], ray[2]);
			float* temp = traceRay(ray, maxDepth);

			unsigned char colToAdd[3] = { 0, 0, 0 };
			colToAdd[0] = (unsigned char)temp[0];
			colToAdd[1] = (unsigned char)temp[1];
			colToAdd[2] = (unsigned char)temp[2];

			out[j][i][0] = colToAdd[0];
			out[j][i][1] = colToAdd[1];
			out[j][i][2] = colToAdd[2];
		}
	}
	float c = 255 / maxCol;
	for (int i = 0; i < outWidth; i++) {
		for (int j = 0; j < outHeight; j++) {
			out[j][i][0] *= c;
			out[j][i][1] *= c;
			out[j][i][2] *= c;
		}
	}

	stbi_write_png("out.png", outWidth, outHeight, channels, out, outWidth * channels * sizeof(unsigned char));
	printf("Wrote to png!");
}