#include <iostream>
#include <cmath>
#include <vector>
#include <string>

struct Coord {
	double x;
	double y;

	/// <summary>
	/// Soma los vectores entre si
	/// </summary>
	/// <param name="b"></param>
	/// <returns></returns>
	Coord operator+ (Coord b) {
		Coord temp = { x + b.x, y + b.y };
		return temp;
	}

	/// <summary>
	/// Se fija si las coordenada de los vectores son las mismas
	/// </summary>
	/// <param name="b"></param>
	/// <returns></returns>
	bool operator == (Coord b) {
		return (x == b.x && y == b.y);
	}

	/// <summary>
	/// Se fija si uno de los valores de las coordenadas coinciden
	/// </summary>
	/// <param name="b"></param>
	/// <returns></returns>
	bool operator |= (Coord b)
	{
		return (x == b.x || y == b.y);
	}
};

struct Vector
{
	int id = 0;
	Coord from;
	Coord to;
};

struct Point {
	int id1;	//primer vector que forma el punto gracias a una interseccion
	int id2;	//segundo vector que forma el punto gracias a una interseccion
	bool isConnected;
	Coord coord;
};

std::vector<Vector> quad;

/// <summary>
/// Analiza si los vectores forman un cuadrilatero
/// </summary>
/// <param name="up"></param>
/// <param name="right"></param>
/// <param name="down"></param>
/// <param name="left"></param>
void IsARectangle(Vector up, Vector right, Vector down, Vector left);

/// <summary>
/// Revisa si los vectores se tocan
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <returns></returns>
Point VectorIntersect(Vector a, Vector b);

/// <summary>
/// Revisa si los vectores son paralelas
/// </summary>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="p3"></param>
/// <param name="p4"></param>
/// <returns></returns>
bool IsParallel(Coord p1, Coord p2, Coord p3, Coord p4);

/// <summary>
/// Revisa si los vectores se tocan en algun punto en el espacio
/// </summary>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="p3"></param>
/// <param name="p4"></param>
/// <returns></returns>
Coord ConnectionPoint(Coord p1, Coord p2, Coord p3, Coord p4);


/// <summary>
/// Encuentra el perimetro sumando la magnitud de sus lados
/// </summary>
/// <param name="leftUp"></param>
/// <param name="rightUp"></param>
/// <param name="leftDown"></param>
/// <param name="rightDown"></param>
/// <returns></returns>
float RectPerimeter(Coord leftUp, Coord rightUp, Coord leftDown, Coord rightDown);

/// <summary>
/// Calcula el area de un cuadrilatero.
/// </summary>
/// <param name="leftUp"></param>
/// <param name="rightUp"></param>
/// <param name="leftDown"></param>
/// <param name="rightDown"></param>
/// <returns></returns>
float RectArea(Coord leftUp, Coord rightUp, Coord leftDown, Coord rightDown);


/// <summary>
/// Ordena los vectores de izquierda a derecha
/// </summary>
/// <param name="vector"></param>
void SortVector(Vector& vector);

void PointComb(std::vector<Point> arr, std::vector<std::vector<Point>>& combResult);

void BuildSidesForEveryCombination(std::vector<std::vector<Point>> combResult);

bool ThisPointsMake4Sides(std::vector<Point> v);

bool ThisSidesMakeARectangle(std::vector<Vector> sides);

void CalculateArea(std::vector<Vector> quad);

float Magnitude(Vector v);

void main()
{

	Vector a = Vector{ 0, Coord{0.0,2.0}, Coord{-0.2f, 1.2f} };
	Vector b = Vector{ 1, Coord{.4f, 1.6f}, Coord{0.2, 0.8f} };

	Vector c = Vector{ 2, Coord{ -6.0, -2.0 }, Coord{ 1.0,3.0 } };
	Vector d = Vector{ 3, Coord{ -4.0,-3.0 }, Coord{ 0.0,3.0 } };



	// Los ordeno para que el principio del vector siempre quede del lado izquierdo y el final del derecho
	SortVector(a);
	SortVector(b);
	SortVector(c);
	SortVector(d);

	IsARectangle(a, b, c, d);

	system("pause");
}

void IsARectangle(Vector up, Vector right, Vector down, Vector left) {


	//Los cuatro vectores utilizados en el ejercicio
	std::vector<Vector> vectors = { up, right, down, left };

	//Los puntos en donde los vectores se superponen
	std::vector<Point> connections;

	//Los vectores resultantes de combinar los puntos
	std::vector<std::vector<Point>> combResult;

	int index = 0;

	for (int i = 0; i < vectors.size(); i++) {	//Recorro todos los vectores dados
		for (int j = index; j < 4; j++) {	//Produzco todas las combinaciones entre estos vectores sin que las combinaciones se repitan
			if (i != j) {	//Mientras que no esté compartando los mismos vectores
				connections.push_back(VectorIntersect(vectors[i], vectors[j]));	//Ingreso en mi vector de puntos, los puntos donde los vectores interesctan (si es que lo hacen).
			}
		}

		index++;
	}

	PointComb(connections, combResult);

	BuildSidesForEveryCombination(combResult);

	CalculateArea(quad);
}

/// <summary>
/// Corrobora si es un cuadrilatero
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <returns></returns>
Point VectorIntersect(Vector a, Vector b) {
	Coord p1 = a.from;
	Coord p2 = a.to;
	Coord p3 = b.from;
	Coord p4 = b.to;

	Point point = Point{ 0, 0, false, Coord{0, 0} };

	// Comprobar si son paralelas
	if (IsParallel(p1, p2, p3, p4)) {
		return point;
	}
	else {
		// Si no son paralelas se busca los determinantes.

		// Busco las determinantes
		point.coord = ConnectionPoint(p1, p2, p3, p4);

		point.id1 = a.id;
		point.id2 = b.id;

		// Me fijo que las coordenadas esten dentro de los vectores originales.
		bool inX1 = (p1.x <= point.coord.x && point.coord.x <= p2.x) || (p1.x >= point.coord.x && point.coord.x >= p2.x);
		bool inX2 = (p3.x <= point.coord.x && point.coord.x <= p4.x) || (p3.x >= point.coord.x && point.coord.x >= p4.x);

		bool inY1 = (p1.y <= point.coord.y && point.coord.y <= p2.y) || (p1.y >= point.coord.y && point.coord.y >= p2.y);
		bool inY2 = (p3.y <= point.coord.y && point.coord.y <= p4.y) || (p3.y >= point.coord.y && point.coord.y >= p4.y);


		// Corroboro que las condiciones se cumplan para saber si es un cuadrilatero.
		if ((inX1 || inX2) && (inY1 || inY2)) {
			point.isConnected = true;
			return point; // Es un cuadrilatero.
		}
		else {
			point.isConnected = false;
			return point; // No es un cuadrilatero.
		}
	}
}

/// <summary>
/// Revisa si son paralelos los vectores
/// </summary>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="p3"></param>
/// <param name="p4"></param>
/// <returns></returns>
bool IsParallel(Coord p1, Coord p2, Coord p3, Coord p4) {	//Cross product. Si el cross product da 0 entonces los vectores son paralelos
	int num = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);
	std::cout << num;
	//system("pause");
	return num == 0;
}


/// <summary>
/// Devuelve la coordenada de conefrom.xon
/// </summary>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="p3"></param>
/// <param name="p4"></param>
/// <returns></returns>
Coord ConnectionPoint(Coord p1, Coord p2, Coord p3, Coord p4) {
	Coord result;


	// punto de interseccion en X
	result.x = ((p1.x * p2.y - p1.y * p2.x) * (p3.x - p4.x) - (p1.x - p2.x) * (p3.x * p4.y - p3.y * p4.x)) / ((p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x));

	// punto de interseccion en Y
	result.y = ((p1.x * p2.y - p1.y * p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x * p4.y - p3.y * p4.x)) / ((p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x));

	return result;
}

/// <summary>
/// Perimetro del cuadrilatero
/// </summary>
/// <param name="leftUp"></param>
/// <param name="rightUp"></param>
/// <param name="leftDown"></param>
/// <param name="rightDown"></param>
/// <returns></returns>
float RectPerimeter(std::vector<Vector> sides) {
	// Busco la magnitud de cada "Cara" para poder sacar el perimetro.
	double side0Magnitude = sqrt(((sides[0].to.x - sides[0].from.x) * (sides[0].to.x - sides[0].from.x)) + ((sides[0].to.y - sides[0].from.y) * (sides[0].to.y - sides[0].from.y)));
	double side1Magnitude = sqrt(((sides[1].to.x - sides[1].from.x) * (sides[1].to.x - sides[1].from.x)) + ((sides[1].to.y - sides[1].from.y) * (sides[1].to.y - sides[1].from.y)));
	double side2Magnitude = sqrt(((sides[2].to.x - sides[2].from.x) * (sides[2].to.x - sides[2].from.x)) + ((sides[2].to.y - sides[2].from.y) * (sides[2].to.y - sides[2].from.y)));
	double side3Magnitude = sqrt(((sides[3].to.x - sides[3].from.x) * (sides[3].to.x - sides[3].from.x)) + ((sides[3].to.y - sides[3].from.y) * (sides[3].to.y - sides[3].from.y)));

	system("cls");

	std::cout << "Lado 0: " << sides[0].from.x << "," << sides[0].from.y << " hasta " << sides[0].to.x << "," << sides[0].to.y << std::endl;
	std::cout << "Lado 1: " << sides[1].from.x << "," << sides[1].from.y << " hasta " << sides[1].to.x << "," << sides[1].to.y << std::endl;
	std::cout << "Lado 2: " << sides[2].from.x << "," << sides[2].from.y << " hasta " << sides[2].to.x << "," << sides[2].to.y << std::endl;
	std::cout << "Lado 3: " << sides[3].from.x << "," << sides[3].from.y << " hasta " << sides[3].to.x << "," << sides[3].to.y << std::endl;

	return side0Magnitude + side1Magnitude + side2Magnitude + side3Magnitude;
}

float RectArea(Coord leftUp, Coord rightUp, Coord leftDown, Coord rightDown) {
//needs fix
	return 0.0f;
}

/// <summary>
/// Ordena las coordenadas del vector de menor a mayor
/// </summary>
/// <param name="vector"></param>
void SortVector(Vector& vector)
{
	if (vector.to.x < vector.from.x) {
		float savex = vector.from.x;
		vector.from.x = vector.to.x;
		vector.to.x = savex;
		float savey = vector.from.y;
		vector.from.y = vector.to.y;
		vector.to.y = savey;
	}
	else if (vector.from.x == vector.to.x) {
		if (vector.to.y < vector.from.y) {
			float savey = vector.from.y;
			vector.from.y = vector.to.y;
			vector.to.y = savey;
		}
	}
}

/// <summary>
/// Calcula todas las posibles combinaciones
/// </summary>
/// <param name="arr"></param>
/// <param name="combResult"></param>
void PointComb(std::vector<Point> arr, std::vector<std::vector<Point>>& combResult) {
	//HACE TODAS LAS FIGURAS POSIBLES DE FORMAR CON LOS PUNTOS QUE TENEMOS DE ConnectionPoint
	for (int i = 0; i < arr.size(); i++) {	//Recorro todo un array con puntos que intersectan y algunos que no intersectan
		for (int j = i + 1; j < arr.size(); j++) {
			std::vector<Point> comb;
			for (int z = 0; z < arr.size(); z++) {
				if (z != j && z != i)
					comb.push_back(arr[z]);
			}

			combResult.push_back(comb);
		}
	}
}

void BuildSidesForEveryCombination(std::vector<std::vector<Point>> combResult) {
	for (int i = 0; i < combResult.size(); i++) {	//Recorro todas las combinaciones
		if (ThisPointsMake4Sides(combResult[i]))	//Chequeo si esta combinacion forma 4 lados
			break;
	}
}

bool ThisPointsMake4Sides(std::vector<Point> v) {
	int z = 0;
	std::vector<Vector> result;
	for (int i = 0; i < v.size(); i++) {	//Recorro los puntos
		for (int j = i + 1; j < v.size(); j++) {
			//Si el vector que formó a X punto es el mismo que formó a Y punto, y ambos puntos son distintos:
			if ((v[i].id1 == v[j].id1 || v[i].id2 == v[j].id2 || v[i].id1 == v[j].id2 || v[i].id2 == v[j].id1) && i != j) {
				Vector side0;
				side0.from = v[i].coord;
				side0.to = v[j].coord;
				side0.id = z;
				z++;
				result.push_back(side0);
			}
		}
		if (result.size() == 4) {
			for (int i = 0; i < 4; i++) {
				SortVector(result[i]);
			}
			if (ThisSidesMakeARectangle(result)) {
				for (int i = 0; i < result.size(); i++) {
					quad.push_back(result[i]);
				}
				std::cout << RectPerimeter(result) << std::endl;
				return true;
			}

		}
	}
	return false;

}

bool ThisSidesMakeARectangle(std::vector<Vector> sides) {
	float totalAngle = 0;
	float angle = 0;
	int count = 0;
	for (int i = 0; i < sides.size(); i++) {
		for (int j = i + 1; j < sides.size(); j++) {
			if (((float)sides[i].from.x == (float)sides[j].from.x && (float)sides[i].from.y == (float)sides[j].from.y) || ((float)sides[i].to.x == (float)sides[j].to.x && (float)sides[i].to.y == (float)sides[j].to.y) || ((float)sides[i].from.x == (float)sides[j].to.x && (float)sides[i].from.y == (float)sides[j].to.y) || ((float)sides[i].to.x == (float)sides[j].from.x && (float)sides[i].to.y == (float)sides[j].from.y)) {
				count++;
			}
		}
	}
	if (count == 4) {
		for (int i = 0; i < sides.size(); i++) {
			for (int j = i + 1; j < sides.size(); j++) {
				if (((float)sides[i].from.x == (float)sides[j].from.x && (float)sides[i].from.y == (float)sides[j].from.y) || ((float)sides[i].to.x == (float)sides[j].to.x && (float)sides[i].to.y == (float)sides[j].to.y) && i != j) {
					double a = (sides[i].to.x - sides[i].from.x) * (sides[j].to.x - sides[j].from.x);
					double b = (sides[i].to.y - sides[i].from.y) * (sides[j].to.y - sides[j].from.y);
					double c = Magnitude(sides[i]);
					double d = Magnitude(sides[j]);
					angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
					if (angle >= 180) {
						break;
					}
					totalAngle += angle;
				}
				else if ((float)sides[i].from.x == (float)sides[j].to.x && (float)sides[i].from.y == (float)sides[j].to.y && i != j) {
					Coord save = sides[j].to;
					sides[j].to = sides[j].from;
					sides[j].from = save;
					double a = (sides[i].to.x - sides[i].from.x) * (sides[j].to.x - sides[j].from.x);
					double b = (sides[i].to.y - sides[i].from.y) * (sides[j].to.y - sides[j].from.y);
					double c = Magnitude(sides[i]);
					double d = Magnitude(sides[j]);
					angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
					totalAngle += angle;
					if (angle >= 180) {
						break;
					}
					save = sides[j].to;
					sides[j].to = sides[j].from;
					sides[j].from = save;
				}
				else if ((float)sides[i].to.x == (float)sides[j].from.x && (float)sides[i].to.y == (float)sides[j].from.y && i != j) {
					Coord save = sides[i].to;
					sides[i].to = sides[i].from;
					sides[i].from = save;
					double a = (sides[i].to.x - sides[i].from.x) * (sides[j].to.x - sides[j].from.x);
					double b = (sides[i].to.y - sides[i].from.y) * (sides[j].to.y - sides[j].from.y);
					double c = Magnitude(sides[i]);
					double d = Magnitude(sides[j]);
					angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
					totalAngle += angle;
					if (angle >= 180) {
						break;
					}
					save = sides[i].to;
					sides[i].to = sides[i].from;
					sides[i].from = save;
				}
			}

		}
		return ((int)totalAngle == 360);
	}
	else {
		return false;
	}

}

void CalculateArea(std::vector<Vector> quad)
{
	Vector newLine;
	for (int j = 1; j < quad.size(); j++) {
		if ((float)quad[0].from.x == (float)quad[j].from.x && (float)quad[0].from.y == (float)quad[j].from.y) {
			newLine.from = quad[0].to;
			newLine.to = quad[j].to;
		}
		else if ((float)quad[0].from.x == (float)quad[j].to.x && (float)quad[0].from.y == (float)quad[j].to.y) {
			newLine.from = quad[0].to;
			newLine.to = quad[j].from;
		}
		else if ((float)quad[0].to.x == (float)quad[j].from.x && (float)quad[0].to.y == (float)quad[j].from.y) {
			newLine.from = quad[0].to;
			newLine.to = quad[j].from;
		}
		else if ((float)quad[0].to.x == (float)quad[j].to.x && (float)quad[0].to.y == (float)quad[j].to.y) {
			newLine.from = quad[0].from;
			newLine.to = quad[j].from;
		}
	}
	std::vector<float> heights;
	float angle = 0;
	float height = 0;
	for (int i = 0; i < quad.size(); i++) {
		if (((float)quad[i].from.x == (float)newLine.from.x && (float)quad[i].from.y == (float)newLine.from.y)) {
			double a = (quad[i].to.x - quad[i].from.x) * (newLine.to.x - newLine.from.x);
			double b = (quad[i].to.y - quad[i].from.y) * (newLine.to.y - newLine.from.y);
			double c = Magnitude(quad[i]);
			double d = Magnitude(newLine);
			angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
			if (angle < 90) {
				height = ((sin(angle)) * c);
				heights.push_back(height);
			}

		}
		else if ((float)quad[i].from.x == (float)newLine.to.x && (float)quad[i].from.y == (float)newLine.to.y) {
			Coord save = newLine.to;
			newLine.to = newLine.from;
			newLine.from = save;
			double a = (quad[i].to.x - quad[i].from.x) * (newLine.to.x - newLine.from.x);
			double b = (quad[i].to.y - quad[i].from.y) * (newLine.to.y - newLine.from.y);
			double c = Magnitude(quad[i]);
			double d = Magnitude(newLine);
			angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
			if (angle < 90) {
				height = ((sin(angle)) * c);
				heights.push_back(height);
			}
			save = newLine.to;
			newLine.to = newLine.from;
			newLine.from = save;
		}
		else if ((float)quad[i].to.x == (float)newLine.from.x && (float)quad[i].to.y == (float)newLine.from.y) {
			Coord save = quad[i].to;
			quad[i].to = quad[i].from;
			quad[i].from = save;
			double a = (quad[i].to.x - quad[i].from.x) * (newLine.to.x - newLine.from.x);
			double b = (quad[i].to.y - quad[i].from.y) * (newLine.to.y - newLine.from.y);
			double c = Magnitude(quad[i]);
			double d = Magnitude(newLine);
			angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
			if (angle < 90) {
				height = ((sin(angle)) * c);
				heights.push_back(height);
			}
			save = quad[i].to;
			quad[i].to = quad[i].from;
			quad[i].from = save;
		}
		else if (((float)quad[i].to.x == (float)newLine.to.x && (float)quad[i].to.y == (float)newLine.to.y)) {
			Coord save = quad[i].to;
			quad[i].to = quad[i].from;
			quad[i].from = save;
			save = newLine.to;
			newLine.to = newLine.from;
			newLine.from = save;
			double a = (quad[i].to.x - quad[i].from.x) * (newLine.to.x - newLine.from.x);
			double b = (quad[i].to.y - quad[i].from.y) * (newLine.to.y - newLine.from.y);
			double c = Magnitude(quad[i]);
			double d = Magnitude(newLine);
			angle = (acos((float)((a + b) / (c * d)))) * 180 / 3.1415;
			if (angle < 90) {
				height = ((sin(angle)) * c);
				heights.push_back(height);
			}
			save = quad[i].to;
			quad[i].to = quad[i].from;
			quad[i].from = save;
			save = newLine.to;
			newLine.to = newLine.from;
			newLine.from = save;
		}
	}
	int i = 0;
	std::cout << "Total area: ";
}


/// <summary>
/// Calcula la magnitud de un vector (el modulo)
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
float Magnitude(Vector v) {
	return sqrt(((v.to.x - v.from.x) * (v.to.x - v.from.x)) + ((v.to.y - v.from.y) * (v.to.y - v.from.y)));
}


