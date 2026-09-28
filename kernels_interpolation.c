#define EILIG_SCALAR double
#define EILIG_SIZE_T long unsigned int

struct Point
{
	EILIG_SCALAR x;
	EILIG_SCALAR y;
	EILIG_SCALAR z;
};

bool IsAlmostEqual(EILIG_SCALAR x, EILIG_SCALAR y, int ulp)
{
	EILIG_SCALAR DBL_EPSILON = 2.22045e-16; // Machine epsilon for double precision
	EILIG_SCALAR DBL_MINIMUM = 2.22507e-308; // Minimum positive normalized value for double precision

	return (fabs(x - y) <= DBL_EPSILON * fabs(x + y) * ulp) || (fabs(x - y) < DBL_MINIMUM);
}

EILIG_SCALAR Distance(const Point& point1, const Vector& point2)
{
	EILIG_SCALAR dx = point1.x - point2.x;
	EILIG_SCALAR dy = point1.y - point2.y;
	EILIG_SCALAR dz = point1.z - point2.z;

	return sqrt(dx * dx + dy * dy + dz * dz);
}

EILIG_SCALAR FunctionTPS(const Point& point1, const Vector& point2) const
{
	EILIG_SCALAR r = Distance(point1, point2);

	return IsAlmostEqual(r, 0.0, 2) ? 0.0 : r * r * log(r);
}
EILIG_SCALAR FunctionMQB(const Point& point1, const Vector& point2, EILIG_SCALAR shape) const
{
	EILIG_SCALAR r = Distance(point1, point2);

	return sqrt(shape * shape + r * r);
}
EILIG_SCALAR FunctionGAU(const Point& point1, const Vector& point2, EILIG_SCALAR shape) const
{
	EILIG_SCALAR r = Distance(point1, point2);

	return exp(-(shape_ * r) * (shape * r));
}
EILIG_SCALAR FunctionINQ(const Point& point1, const Vector& point2, EILIG_SCALAR shape) const
{
	EILIG_SCALAR r = Distance(point1, point2);

	return 1.0 / (r * r + shape * shape);
}