#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"
#include "performance_test.h"

// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
	// 为当前帧启动性能指标收集
	renderer.startFrameMetrics();

#if defined(USE_SIMD_OPTIMIZATION)
	// SIMD优化版本
	MatrixSIMD p = renderer.perspective * MatrixSIMD(camera) * MatrixSIMD(mesh->world);

	for (triIndices& ind : mesh->triangles) {
		renderer.recordTriangleProcessed();  // 记录三角形处理

		VertexSIMD t[3];
		for (unsigned int i = 0; i < 3; i++) {
			t[i].p = p * Vec4SIMD(mesh->vertices[ind.v[i]].p);
			t[i].p.divideW();

			t[i].normal = MatrixSIMD(mesh->world) * Vec4SIMD(mesh->vertices[ind.v[i]].normal);
			t[i].normal.normalise();

			// 转换到屏幕空间
			renderer.applyViewportTransform(t[i].p);

			t[i].rgb = mesh->vertices[ind.v[i]].rgb;
		}

		// 检查Z值是否在有效范围内
		if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f)
			continue;

		triangle tri(t[0], t[1], t[2]);
		tri.draw(renderer, L, mesh->ka, mesh->kd);
	}

#else 
	// Combine perspective, camera, and world transformations for the mesh
	matrix p = renderer.perspective * camera * mesh->world;

	// Iterate through all triangles in the mesh
	for (triIndices& ind : mesh->triangles) {
		Vertex t[3]; // Temporary array to store transformed triangle vertices

		// Transform each vertex of the triangle
		for (unsigned int i = 0; i < 3; i++) {
			t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
			t[i].p.divideW(); // Perspective division to normalize coordinates

			// Transform normals into world space for accurate lighting
			// no need for perspective correction as no shearing or non-uniform scaling
			t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
			t[i].normal.normalise();

			// Map normalized device coordinates to screen space
			t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
			t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
			t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

			// Copy vertex colours
			t[i].rgb = mesh->vertices[ind.v[i]].rgb;
		}

		// Clip triangles with Z-values outside [-1, 1]
		if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

		// Create a triangle object and render it
		triangle tri(t[0], t[1], t[2]);
		tri.draw(renderer, L, mesh->ka, mesh->kd);
	}
#endif
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
	Renderer renderer;
	// create light source {direction, diffuse intensity, ambient intensity}
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
	// camera is just a matrix
	matrix camera = matrix::makeIdentity(); // Initialize the camera with identity matrix

	bool running = true; // Main loop control variable

	std::vector<Mesh*> scene; // Vector to store scene objects

	// Create a sphere and a rectangle mesh
	Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
	//Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

	// add meshes to scene
	scene.push_back(&mesh);
	// scene.push_back(&mesh2); 

	float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
	mesh.world = matrix::makeTranslation(x, y, z);
	//mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);

	// Main rendering loop
	while (running) {
		renderer.canvas.checkInput(); // Handle user input
		renderer.clear(); // Clear the canvas for the next frame

		// Apply transformations to the meshes
	 //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
		mesh.world = matrix::makeTranslation(x, y, z);

		// Handle user inputs for transformations
		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
		if (renderer.canvas.keyPressed('A')) x += -0.1f;
		if (renderer.canvas.keyPressed('D')) x += 0.1f;
		if (renderer.canvas.keyPressed('W')) y += 0.1f;
		if (renderer.canvas.keyPressed('S')) y += -0.1f;
		if (renderer.canvas.keyPressed('Q')) z += 0.1f;
		if (renderer.canvas.keyPressed('E')) z += -0.1f;

		// Render each object in the scene
		for (auto& m : scene)
			render(renderer, m, camera, L);

		renderer.present(); // Display the rendered frame
	}
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
	unsigned int r = rng.getRandomInt(0, 3);

	switch (r) {
	case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
	case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
	case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
	default: return matrix::makeIdentity();
	}
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
	Renderer renderer;
	matrix camera;
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	bool running = true;

	std::vector<Mesh*> scene;

	// Create a scene of 40 cubes with random rotations
	for (unsigned int i = 0; i < 20; i++) {
		Mesh* m = new Mesh();
		*m = Mesh::makeCube(1.f);
		m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
		scene.push_back(m);
		m = new Mesh();
		*m = Mesh::makeCube(1.f);
		m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
		scene.push_back(m);
	}

	float zoffset = 8.0f; // Initial camera Z-offset
	float step = -0.1f;  // Step size for camera movement

	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;

	// Main rendering loop
	while (running) {
		renderer.canvas.checkInput();
		renderer.clear();

		camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

		// Rotate the first two cubes in the scene
		scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
		scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		zoffset += step;
		if (zoffset < -60.f || zoffset > 8.f) {
			step *= -1.f;
			if (++cycle % 2 == 0) {
				end = std::chrono::high_resolution_clock::now();
				std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
				start = std::chrono::high_resolution_clock::now();
			}
		}

		for (auto& m : scene)
			render(renderer, m, camera, L);
		renderer.present();
	}

	for (auto& m : scene)
		delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
	Renderer renderer;
	matrix camera = matrix::makeIdentity();
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	std::vector<Mesh*> scene;

	struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
	std::vector<rRot> rotations;

	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

	// Create a grid of cubes with random rotations
	for (unsigned int y = 0; y < 6; y++) {
		for (unsigned int x = 0; x < 8; x++) {
			Mesh* m = new Mesh();
			*m = Mesh::makeCube(1.f);
			scene.push_back(m);
			m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
			rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
			rotations.push_back(r);
		}
	}

	// Create a sphere and add it to the scene
	Mesh* sphere = new Mesh();
	*sphere = Mesh::makeSphere(1.0f, 10, 20);
	scene.push_back(sphere);
	float sphereOffset = -6.f;
	float sphereStep = 0.1f;
	sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;

	bool running = true;
	while (running) {
		renderer.canvas.checkInput();
		renderer.clear();

		// Rotate each cube in the grid
		for (unsigned int i = 0; i < rotations.size(); i++)
			scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

		// Move the sphere back and forth
		sphereOffset += sphereStep;
		sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
		if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
			sphereStep *= -1.f;
			if (++cycle % 2 == 0) {
				end = std::chrono::high_resolution_clock::now();
				std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
				start = std::chrono::high_resolution_clock::now();
			}
		}

		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		for (auto& m : scene)
			render(renderer, m, camera, L);
		renderer.present();
	}

	for (auto& m : scene)
		delete m;
}

void scene3() {
	Renderer renderer;
	matrix camera = matrix::makeIdentity();
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	std::vector<Mesh*> scene;

	// 设置网格参数
	const int GRID_SIZE = 8;     // 8x8的网格提供了64个球体，这个数量足够产生明显的计算负载
	const float SPACING = 2.0f;   // 球体间距设置为2.0保证有适量的重叠
	const int SPHERE_DETAIL = 20; // 球体的细分程度，影响每个球体的顶点数量

	// 创建球体网格
	for (int y = 0; y < GRID_SIZE; y++) {
		for (int x = 0; x < GRID_SIZE; x++) {
			Mesh* sphere = new Mesh();
			// 使用较高的细分度创建球体，增加计算负载
			*sphere = Mesh::makeSphere(1.0f, SPHERE_DETAIL, SPHERE_DETAIL * 2);

			// 计算球体在网格中的位置
			float xPos = (x - GRID_SIZE / 2) * SPACING;
			float yPos = (y - GRID_SIZE / 2) * SPACING;
			sphere->world = matrix::makeTranslation(xPos, yPos, -10.0f);

			// 给每个球体设置稍微不同的颜色，便于视觉区分
			float r = 0.5f + (float)x / GRID_SIZE * 0.5f;
			float g = 0.5f + (float)y / GRID_SIZE * 0.5f;
			float b = 0.7f;
			sphere->setColour(colour(r, g, b), 0.2f, 0.8f);

			scene.push_back(sphere);
		}
	}

	/*auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;*/

	int frameCount = 0;
	auto totalStart = std::chrono::high_resolution_clock::now();

	bool running = true;
	while (running) {
		renderer.canvas.checkInput();
		renderer.clear();

		// 为每个球体创建独特的动画
		for (size_t i = 0; i < scene.size(); i++) {
			// 使用球体在网格中的位置计算旋转角度，创造波浪效果
			int gridX = i % GRID_SIZE;
			int gridY = i / GRID_SIZE;
			float timeOffset = static_cast<float>(frameCount) * 0.01f;
			float xAngle = std::sin(timeOffset + gridX * 0.5f) * 0.02f;
			float yAngle = std::cos(timeOffset + gridY * 0.5f) * 0.02f;

			scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(xAngle, yAngle, 0.0f);
		}

		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		// 渲染场景
		for (auto& m : scene)
			render(renderer, m, camera, L);

		renderer.present();
		frameCount++;

		// 每100帧输出一次性能报告
		if (frameCount % 100 == 0) {
			std::cout << "Frame " << frameCount << "\n";
			std::cout << renderer.getPerformanceReport();
			std::cout << "-------------------\n";
		}
	}
	// 输出总体性能统计
	auto totalEnd = std::chrono::high_resolution_clock::now();
	double totalTime = std::chrono::duration<double, std::milli>(totalEnd - totalStart).count();

	std::cout << "\nFinal Performance Statistics:\n";
	std::cout << "Total Frames: " << frameCount << "\n";
	std::cout << "Total Time: " << totalTime << "ms\n";
	std::cout << "Average FPS: " << (frameCount * 1000.0 / totalTime) << "\n";

	// 清理资源
	for (auto& m : scene)
		delete m;

}

// Entry point of the application
// No input variables
int main() {
	// Uncomment the desired scene function to run
	//scene1();
	//scene2();
	//sceneTest(); 


	PerformanceBenchmark benchmark;

	while (true) {
		std::cout << "\n=== Rasterization Testing System 1.0 ===\n";
		std::cout << "1. Scene1 test (cubes)\n";
		std::cout << "2. Scene2 test (cube mesh)\n";
		std::cout << "3. Scene3 test (sphere mesh)\n";
		std::cout << "4. simple sence test\n";
		std::cout << "5. RUN THE FULL PERFORMANCE TEST!(All SCENE)\n";
		std::cout << "0. Quit\n";

		int choice;
		std::cin >> choice;

		switch (choice) {
		case 0:
			return 0;
		case 1:
		case 2:
		case 3:
		case 4: {
			// 手动测试单个场景
			void (*sceneFunc)() = nullptr;
			switch (choice) {
			case 1: sceneFunc = scene1; break;
			case 2: sceneFunc = scene2; break;
			case 3: sceneFunc = scene3; break;
			case 4: sceneFunc = sceneTest; break;
			}
			sceneFunc();
			break;
		}
		case 5: {
			// 运行完整性能测试
			std::cout << "Start performance test...\n";

			// Scene1测试
			std::cout << "Scene1...\n";
#define USE_SIMD_OPTIMIZATION 0
			benchmark.runBenchmark("Scene1", scene1, false);
#undef USE_SIMD_OPTIMIZATION
#define USE_SIMD_OPTIMIZATION 1
			benchmark.runBenchmark("Scene1", scene1, true);
#undef USE_SIMD_OPTIMIZATION

			// Scene2测试
			std::cout << "Scene2...\n";
#define USE_SIMD_OPTIMIZATION 0
			benchmark.runBenchmark("Scene2", scene2, false);
#undef USE_SIMD_OPTIMIZATION
#define USE_SIMD_OPTIMIZATION 1
			benchmark.runBenchmark("Scene2", scene2, true);
#undef USE_SIMD_OPTIMIZATION

			// Scene3测试
			std::cout << "Scene3...\n";
#define USE_SIMD_OPTIMIZATION 0
			benchmark.runBenchmark("Scene3", scene3, false);
#undef USE_SIMD_OPTIMIZATION
#define USE_SIMD_OPTIMIZATION 1
			benchmark.runBenchmark("Scene3", scene3, true);
#undef USE_SIMD_OPTIMIZATION

			// 生成报告
			benchmark.generateReport("performance_report.txt");
			std::cout << "Performance Test Finished! the report has been written in performance_report.txt\n";
			break;
		}
		default:
			std::cout << "Invalid Choice!\n";
		}
	}
	return 0;
}