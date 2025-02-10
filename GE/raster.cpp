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


//// 定义优化模式枚举
//enum OptimizationMode {
//	BASE,           // 基础版本
//	SIMD,           // SIMD优化版本
//	THREADED,       // 多线程优化版本
//	SIMD_THREADED   // SIMD+多线程组合优化
//};

int g_maxCycles;
bool g_isPerformanceTest;

// 设置场景运行模式的函数
void setSceneRunMode(bool isPerformanceTest) {
	if (isPerformanceTest) {
		// 性能测试模式：限制cycle数
		g_maxCycles = 1;  // 可以根据需要调整这个值
		g_isPerformanceTest = true;
	}
	else {
		// 单一场景模式：无限制运行
		g_maxCycles = -1;  // 使用-1表示无限制
		g_isPerformanceTest = false;
	}
}

// 主要渲染函数，支持不同的优化模式
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L, OptimizationType mode) {
	// 启动性能指标收集
	renderer.startFrameMetrics();

	// 组合所有变换矩阵
	matrix p = renderer.perspective * camera * mesh->world;

	// 处理网格中的每个三角形
	for (triIndices& ind : mesh->triangles) {
		renderer.recordTriangleProcessed();

		Vertex t[3];
		// 变换每个顶点
		for (unsigned int i = 0; i < 3; i++) {
			t[i].p = p * mesh->vertices[ind.v[i]].p;
			t[i].p.divideW();

			t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
			t[i].normal.normalise();

			// 应用视口变换
			renderer.applyViewportTransform(t[i].p);
			t[i].rgb = mesh->vertices[ind.v[i]].rgb;
		}

		// 深度裁剪
		if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f)
			continue;

		// 根据优化模式选择渲染方式
		triangle tri(t[0], t[1], t[2]);
		switch (mode) {
		case BASE:
			tri.draw(renderer, L, mesh->ka, mesh->kd);
			break;
		case SIMD:
#if defined(USE_SIMD_OPTIMIZATION)
			tri.draw(renderer, L, mesh->ka, mesh->kd);
#else
			tri.draw(renderer, L, mesh->ka, mesh->kd);
#endif
			break;
		case THREADED:
			tri.drawThreaded(renderer, L, mesh->ka, mesh->kd);
			break;
		case SIMD_THREADED:
#if defined(USE_SIMD_OPTIMIZATION)
			tri.drawThreaded(renderer, L, mesh->ka, mesh->kd);
#else
			tri.drawThreaded(renderer, L, mesh->ka, mesh->kd);
#endif
			break;
		}
	}
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
			render(renderer, m, camera, L,BASE);

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

	// 性能统计变量初始化 - 这些变量将用于跟踪整个运行过程的性能
	int totalTrianglesProcessed = 0;
	int totalPixelsProcessed = 0;
	double totalFrameTime = 0.0;
	int frameCount = 0;

	// 场景状态控制
	bool running = true;
	std::vector<Mesh*> scene;

	// 创建40个立方体的场景
	for (unsigned int i = 0; i < 20; i++) {
		// 创建并放置左侧立方体
		Mesh* m = new Mesh();
		*m = Mesh::makeCube(1.f);
		m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
		scene.push_back(m);

		// 创建并放置右侧立方体
		m = new Mesh();
		*m = Mesh::makeCube(1.f);
		m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
		scene.push_back(m);
	}

	// 相机和动画控制参数
	float zoffset = 8.0f;
	float step = -0.1f;

	// 时间和周期计数初始化
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;        // 用于计算来回运动
	int cycleNumber = 0;  // 用于显示周期编号

	// 主渲染循环
	while (running) {
		// 记录帧开始时间，用于计算单帧渲染时间
		auto frameStart = std::chrono::high_resolution_clock::now();

		// 基础渲染准备
		renderer.canvas.checkInput();
		renderer.clear();

		// 更新相机位置
		camera = matrix::makeTranslation(0, 0, -zoffset);

		// 更新场景中前两个立方体的旋转
		scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
		scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

		// 检查ESC退出，如果退出也要输出完整统计信息
		if (renderer.canvas.keyPressed(VK_ESCAPE)) {
			std::cout << "\nPerformance Summary for Scene1 (ESC exit):\n";
			std::cout << "Total Cycles Completed: " << cycleNumber << "\n";
			std::cout << "Total Frames: " << frameCount << "\n";
			std::cout << "Average Frame Time: " << totalFrameTime / frameCount << "ms\n";
			std::cout << "Average FPS: " << 1000.0 / (totalFrameTime / frameCount) << "\n";
			std::cout << "Total Triangles Processed: " << totalTrianglesProcessed << "\n";
			std::cout << "Average Triangles per Frame: " << totalTrianglesProcessed / frameCount << "\n";
			std::cout << "Total Pixels Processed: " << totalPixelsProcessed << "\n";
			std::cout << "Average Pixels per Frame: " << totalPixelsProcessed / frameCount << "\n";
			break;
		}

		// 渲染所有物体
		for (auto& m : scene) {
			render(renderer, m, camera, L,BASE);
		}
		renderer.present();

		// 收集每帧的性能数据
		auto frameEnd = std::chrono::high_resolution_clock::now();
		double frameTime = std::chrono::duration<double, std::milli>(frameEnd - frameStart).count();
		totalFrameTime += frameTime;
		frameCount++;
		totalTrianglesProcessed += renderer.getPerformanceMetrics().triangleCount;
		totalPixelsProcessed += renderer.getPerformanceMetrics().pixelsProcessed;

		// 更新相机位置并检查周期完成情况
		zoffset += step;
		if (zoffset < -60.f || zoffset > 8.f) {
			step *= -1.f;
			if (++cycle % 2 == 0) {  // 每完成一个来回运动
				cycleNumber++;        // 增加周期计数
				end = std::chrono::high_resolution_clock::now();

				// 输出当前周期的时间和阶段性统计
				std::cout << cycleNumber << " :"
					<< std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";

				// 在性能测试模式下检查是否达到最大周期数
				if (g_isPerformanceTest && cycleNumber >= g_maxCycles) {
					// 输出最终性能统计
					std::cout << "\nPerformance Summary for Scene1 (Performance Test):\n";
					std::cout << "Total Cycles Completed: " << cycleNumber << "\n";
					std::cout << "Total Frames: " << frameCount << "\n";
					std::cout << "Average Frame Time: " << totalFrameTime / frameCount << "ms\n";
					std::cout << "Average FPS: " << 1000.0 / (totalFrameTime / frameCount) << "\n";
					std::cout << "Total Triangles Processed: " << totalTrianglesProcessed << "\n";
					std::cout << "Average Triangles per Frame: " << totalTrianglesProcessed / frameCount << "\n";
					std::cout << "Total Pixels Processed: " << totalPixelsProcessed << "\n";
					std::cout << "Average Pixels per Frame: " << totalPixelsProcessed / frameCount << "\n";
					running = false;
					break;
				}

				start = std::chrono::high_resolution_clock::now();
			}
		}
	}

	// 清理场景资源
	for (auto& m : scene) {
		delete m;
	}
}
// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
	Renderer renderer;
	matrix camera = matrix::makeIdentity();
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	// 性能统计变量初始化
	int totalTrianglesProcessed = 0;
	int totalPixelsProcessed = 0;
	double totalFrameTime = 0.0;
	int frameCount = 0;

	std::vector<Mesh*> scene;
	struct rRot { float x; float y; float z; };
	std::vector<rRot> rotations;

	// 初始化随机数生成器
	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

	// 创建6x8的立方体网格
	for (unsigned int y = 0; y < 6; y++) {
		for (unsigned int x = 0; x < 8; x++) {
			Mesh* m = new Mesh();
			*m = Mesh::makeCube(1.f);
			scene.push_back(m);
			m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f),
				5.0f - (static_cast<float>(y) * 2.f), -8.f);
			// 为每个立方体生成随机旋转参数
			rRot r{ rng.getRandomFloat(-.1f, .1f),
				   rng.getRandomFloat(-.1f, .1f),
				   rng.getRandomFloat(-.1f, .1f) };
			rotations.push_back(r);
		}
	}

	// 创建并添加球体
	Mesh* sphere = new Mesh();
	*sphere = Mesh::makeSphere(1.0f, 10, 20);
	scene.push_back(sphere);
	float sphereOffset = -6.f;
	float sphereStep = 0.1f;
	sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

	// 时间和周期计数初始化
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;
	int cycleNumber = 0;

	bool running = true;
	while (running) {
		// 记录帧开始时间
		auto frameStart = std::chrono::high_resolution_clock::now();

		renderer.canvas.checkInput();
		renderer.clear();

		// 更新立方体网格的旋转
		for (unsigned int i = 0; i < rotations.size(); i++) {
			scene[i]->world = scene[i]->world *
				matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
		}

		// 更新球体位置
		sphereOffset += sphereStep;
		sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

		// 检查ESC退出并输出统计信息
		if (renderer.canvas.keyPressed(VK_ESCAPE)) {
			std::cout << "\nPerformance Summary for Scene2 (ESC exit):\n";
			std::cout << "Total Cycles Completed: " << cycleNumber << "\n";
			std::cout << "Total Frames: " << frameCount << "\n";
			std::cout << "Average Frame Time: " << totalFrameTime / frameCount << "ms\n";
			std::cout << "Average FPS: " << 1000.0 / (totalFrameTime / frameCount) << "\n";
			std::cout << "Total Triangles Processed: " << totalTrianglesProcessed << "\n";
			std::cout << "Average Triangles per Frame: " << totalTrianglesProcessed / frameCount << "\n";
			std::cout << "Total Pixels Processed: " << totalPixelsProcessed << "\n";
			std::cout << "Average Pixels per Frame: " << totalPixelsProcessed / frameCount << "\n";
			break;
		}

		// 渲染场景中的所有物体
		for (auto& m : scene) {
			render(renderer, m, camera, L,BASE);
		}
		renderer.present();

		// 收集每帧的性能数据
		auto frameEnd = std::chrono::high_resolution_clock::now();
		double frameTime = std::chrono::duration<double, std::milli>(frameEnd - frameStart).count();
		totalFrameTime += frameTime;
		frameCount++;
		totalTrianglesProcessed += renderer.getPerformanceMetrics().triangleCount;
		totalPixelsProcessed += renderer.getPerformanceMetrics().pixelsProcessed;

		// 检查球体运动周期和性能统计
		if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
			sphereStep *= -1.f;
			if (++cycle % 2 == 0) {
				cycleNumber++;
				end = std::chrono::high_resolution_clock::now();

				// 输出当前周期的时间
				std::cout << cycleNumber << " :"
					<< std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";

				// 在性能测试模式下检查是否达到最大周期数
				if (g_isPerformanceTest && cycleNumber >= g_maxCycles) {
					// 输出最终性能统计
					std::cout << "\nPerformance Summary for Scene2 (Performance Test):\n";
					std::cout << "Total Cycles Completed: " << cycleNumber << "\n";
					std::cout << "Total Frames: " << frameCount << "\n";
					std::cout << "Average Frame Time: " << totalFrameTime / frameCount << "ms\n";
					std::cout << "Average FPS: " << 1000.0 / (totalFrameTime / frameCount) << "\n";
					std::cout << "Total Triangles Processed: " << totalTrianglesProcessed << "\n";
					std::cout << "Average Triangles per Frame: " << totalTrianglesProcessed / frameCount << "\n";
					std::cout << "Total Pixels Processed: " << totalPixelsProcessed << "\n";
					std::cout << "Average Pixels per Frame: " << totalPixelsProcessed / frameCount << "\n";
					running = false;
					break;
				}

				start = std::chrono::high_resolution_clock::now();
			}
		}
	}

	// 清理场景资源
	for (auto& m : scene) {
		delete m;
	}
}

void scene3() {
	Renderer renderer;
	matrix camera = matrix::makeIdentity();
	// 调整光照方向和强度，使球体渲染更自然
	Light L{ vec4(0.5f, 0.5f, 1.0f, 0.0f), colour(1.0f, 1.0f, 1.0f), colour(0.4f, 0.4f, 0.4f) };

	// 性能统计变量初始化
	int totalTrianglesProcessed = 0;
	int totalPixelsProcessed = 0;
	double totalFrameTime = 0.0;
	int frameCount = 0;

	std::vector<Mesh*> scene;

	// 创建原子核 - 增加细分度以获得更平滑的效果
	Mesh* nucleus = new Mesh();
	*nucleus = Mesh::makeSphere(1.0f, 64, 128);  // 增加细分度
	nucleus->col.set(1.0f, 0.2f, 0.2f);        // 明亮的红色
	nucleus->ka = 0.6f;                         // 增加环境光系数
	nucleus->kd = 0.8f;                         // 适度的漫反射
	nucleus->world = matrix::makeTranslation(0.0f, 0.0f, -10.0f);
	scene.push_back(nucleus);

	// 创建第一层电子（2个）
	for (int i = 0; i < 2; i++) {
		Mesh* electron = new Mesh();
		*electron = Mesh::makeSphere(0.25f, 16, 32);  // 稍微减小电子大小，增加细分
		electron->col.set(0.3f, 0.7f, 1.0f);         // 明亮的蓝色
		electron->ka = 0.6f;
		electron->kd = 0.8f;
		scene.push_back(electron);
	}

	// 创建第二层电子（6个）
	for (int i = 0; i < 6; i++) {
		Mesh* electron = new Mesh();
		*electron = Mesh::makeSphere(0.25f, 16, 32);
		electron->col.set(0.3f, 1.0f, 0.3f);         // 明亮的绿色
		electron->ka = 0.6f;
		electron->kd = 0.8f;
		scene.push_back(electron);
	}

	// 创建轨道可视化点
	const int ORBIT_POINTS = 48;  // 增加轨道点数量使轨道更平滑
	// 内层轨道
	for (int i = 0; i < ORBIT_POINTS; i++) {
		Mesh* point = new Mesh();
		*point = Mesh::makeSphere(0.04f, 8, 16);     // 更小的轨道点
		point->col.set(0.4f, 0.6f, 0.8f);           // 柔和的蓝色
		point->ka = 0.4f;                           // 减小轨道点的亮度
		point->kd = 0.6f;
		scene.push_back(point);
	}

	// 外层轨道
	for (int i = 0; i < ORBIT_POINTS; i++) {
		Mesh* point = new Mesh();
		*point = Mesh::makeSphere(0.04f, 8, 16);
		point->col.set(0.4f, 0.8f, 0.4f);           // 柔和的绿色
		point->ka = 0.4f;
		point->kd = 0.6f;
		scene.push_back(point);
	}

	// 轨道参数
	const float INNER_RADIUS = 2.8f;  
	const float OUTER_RADIUS = 4.8f;

	// 动画控制
	float animationTime = 0.0f;
	const float CYCLE_DURATION = 8.0f * M_PI; 
	int cycleNumber = 0;

	// 电子速度变化因子
	std::vector<float> electronSpeeds;
	for (int i = 0; i < 8; i++) {
		// 降低速度系数，使动画变慢
		electronSpeeds.push_back(0.2f + RandomNumberGenerator::getInstance().getRandomFloat(-0.2f, 0.2f));
	}

	// 时间计数初始化
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;

	bool running = true;
	while (running) {
		auto frameStart = std::chrono::high_resolution_clock::now();

		renderer.canvas.checkInput();
		renderer.clear();

		// 原子核旋转（降低旋转速度）
		scene[0]->world = matrix::makeTranslation(0.0f, 0.0f, -10.0f) *
			matrix::makeRotateY(animationTime * 0.12f);

		// 更新轨道点位置
		int orbitStartIndex = 8;

		// 内层轨道点
		for (int i = 0; i < ORBIT_POINTS; i++) {
			float angle = (2.0f * M_PI * i) / ORBIT_POINTS;
			float x = INNER_RADIUS * cos(angle);
			float y = INNER_RADIUS * sin(angle);
			scene[orbitStartIndex + i]->world = matrix::makeTranslation(x, y, -10.0f);
		}

		// 外层轨道点
		for (int i = 0; i < ORBIT_POINTS; i++) {
			float angle = (2.0f * M_PI * i) / ORBIT_POINTS;
			float x = OUTER_RADIUS * cos(angle);
			float y = OUTER_RADIUS * sin(angle);
			scene[orbitStartIndex + ORBIT_POINTS + i]->world =
				matrix::makeTranslation(x, y, -10.0f);
		}

		// 更新内层电子位置
		for (int i = 0; i < 2; i++) {
			float angle = animationTime * electronSpeeds[i] + (M_PI * i);
			float x = INNER_RADIUS * cos(angle);
			float y = INNER_RADIUS * sin(angle);
			scene[1 + i]->world = matrix::makeTranslation(x, y, -10.0f);
		}

		// 更新外层电子位置
		for (int i = 0; i < 6; i++) {
			float angle = animationTime * electronSpeeds[i + 2] + (2.0f * M_PI * i / 6.0f);
			float x = OUTER_RADIUS * cos(angle);
			float y = OUTER_RADIUS * sin(angle);
			scene[3 + i]->world = matrix::makeTranslation(x, y, -10.0f);
		}

		// ESC退出检查
		if (renderer.canvas.keyPressed(VK_ESCAPE)) {
			std::cout << "\nPerformance Summary for Scene3 (ESC exit):\n";
			std::cout << "Total Cycles Completed: " << cycleNumber << "\n";
			std::cout << "Total Frames: " << frameCount << "\n";
			std::cout << "Average Frame Time: " << totalFrameTime / frameCount << "ms\n";
			std::cout << "Average FPS: " << 1000.0 / (totalFrameTime / frameCount) << "\n";
			std::cout << "Total Triangles Processed: " << totalTrianglesProcessed << "\n";
			std::cout << "Average Triangles per Frame: " << totalTrianglesProcessed / frameCount << "\n";
			std::cout << "Total Pixels Processed: " << totalPixelsProcessed << "\n";
			std::cout << "Average Pixels per Frame: " << totalPixelsProcessed / frameCount << "\n";
			break;
		}

		// 渲染场景
		for (auto& m : scene) {
			render(renderer, m, camera, L,BASE);
		}
		renderer.present();

		// 性能数据收集
		auto frameEnd = std::chrono::high_resolution_clock::now();
		double frameTime = std::chrono::duration<double, std::milli>(frameEnd - frameStart).count();
		totalFrameTime += frameTime;
		frameCount++;
		totalTrianglesProcessed += renderer.getPerformanceMetrics().triangleCount;
		totalPixelsProcessed += renderer.getPerformanceMetrics().pixelsProcessed;

		// 更新动画时间（可以根据需要调整）
		animationTime += 0.02f;  // 降低时间增量，使动画更慢
		if (animationTime >= CYCLE_DURATION) {
			cycleNumber++;
			end = std::chrono::high_resolution_clock::now();

			std::cout << cycleNumber << " :"
				<< std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";

			if (g_isPerformanceTest && cycleNumber >= g_maxCycles) {
				std::cout << "\nPerformance Summary for Scene3 (Performance Test):\n";
				std::cout << "Total Cycles Completed: " << cycleNumber << "\n";
				std::cout << "Total Frames: " << frameCount << "\n";
				std::cout << "Average Frame Time: " << totalFrameTime / frameCount << "ms\n";
				std::cout << "Average FPS: " << 1000.0 / (totalFrameTime / frameCount) << "\n";
				std::cout << "Total Triangles Processed: " << totalTrianglesProcessed << "\n";
				std::cout << "Average Triangles per Frame: " << totalTrianglesProcessed / frameCount << "\n";
				std::cout << "Total Pixels Processed: " << totalPixelsProcessed << "\n";
				std::cout << "Average Pixels per Frame: " << totalPixelsProcessed / frameCount << "\n";
				running = false;
				break;
			}

			animationTime = 0.0f;
			start = std::chrono::high_resolution_clock::now();
		}
	}

	// 清理资源
	for (auto& m : scene) {
		delete m;
	}
}

// Entry point of the application
// No input variables
int main() {
	PerformanceBenchmark benchmark;

	while (true) {
		std::cout << "\n=== Rasterization Testing System 3.0 ===\n";
		std::cout << "1. Scene1 test (moving cubes)\n";
		std::cout << "2. Scene2 test (cube grid)\n";
		std::cout << "3. Scene3 test (atomic model)\n";
		std::cout << "4. Simple scene test\n";
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
			// 单一场景测试：设置为无限制模式
			setSceneRunMode(false);

			// 选择优化模式
			std::cout << "\nSelect optimization mode:\n";
			std::cout << "1. Base Version\n";
			std::cout << "2. SIMD Optimization\n";
			std::cout << "3. Multi-threaded\n";
			std::cout << "4. SIMD + Multi-threaded\n";

			int modeChoice;
			std::cin >> modeChoice;

			OptimizationType mode;
			switch (modeChoice) {
			case 1: mode = OptimizationType::BASE; break;
			case 2: mode = OptimizationType::SIMD; break;
			case 3: mode = OptimizationType::THREADED; break;
			case 4: mode = OptimizationType::SIMD_THREADED; break;
			default: continue;
			}

			// 根据场景选择运行相应的函数
			switch (choice) {
			case 1:
				scene1();
				break;
			case 2:
				scene2();
				break;
			case 3:
				scene3();
				break;
			case 4:
				sceneTest();
				break;
			}
			break;
		}

		case 5: {

			// 运行完整性能测试
			setSceneRunMode(true);
			
			std::cout << "Start performance test...\n";

			// Scene1测试
			std::cout << "\n\nScene1 Base version...\n";
			benchmark.runBenchmark("Scene1", scene1, OptimizationType::BASE);

			std::cout << "\nScene1 SIMD version...\n";
			benchmark.runBenchmark("Scene1", scene1, OptimizationType::SIMD);

			std::cout << "\nScene1 Threaded version...\n";
			benchmark.runBenchmark("Scene1", scene1, OptimizationType::THREADED);

			std::cout << "\nScene1 SIMD+Thread version...\n";
			benchmark.runBenchmark("Scene1", scene1, OptimizationType::SIMD_THREADED);

			// Scene2测试（类似的四个版本）...
			 std::cout << "\n\nScene2 Base version...\n";
			benchmark.runBenchmark("Scene2", scene2, OptimizationType::BASE);

			std::cout << "\nScene2 SIMD version...\n";
			benchmark.runBenchmark("Scene2", scene2, OptimizationType::SIMD);

			std::cout << "\nScene2 Threaded version...\n";
			benchmark.runBenchmark("Scene2", scene2, OptimizationType::THREADED);

			std::cout << "\nScene3 SIMD+Thread version...\n";
			benchmark.runBenchmark("Scene2", scene2, OptimizationType::SIMD_THREADED);

			// Scene3测试（类似的四个版本）...

			std::cout << "\n\nScene3 Base version...\n";
			benchmark.runBenchmark("Scene3", scene3, OptimizationType::BASE);

			std::cout << "\nScene3 SIMD version...\n";
			benchmark.runBenchmark("Scene3", scene3, OptimizationType::SIMD);

			std::cout << "\nScene3 Threaded version...\n";
			benchmark.runBenchmark("Scene3", scene3, OptimizationType::THREADED);

			std::cout << "Scene3 SIMD+Thread version...\n";
			benchmark.runBenchmark("Scene3", scene3, OptimizationType::SIMD_THREADED);


			// 生成完整报告
			benchmark.generateReport("performance_report.txt");
			std::cout << "Performance Test Finished! Report written to performance_report.txt\n";
			break;
		}

		default:
			std::cout << "Invalid Choice!\n";
		}
	}
	return 0;
}