#include "./includes.hpp"
#include "./utils.hpp"

std::vector<double> interpolationValues;

std::vector<double> getEquidistantPoints(
    double xMin,
    double xMax,
    int valuesCount
) {
    std::vector<double> points(valuesCount);

    double h = (xMax - xMin) / (valuesCount - 1);

    for (int i = 0; i < valuesCount; ++i)
        points[i] = xMin + i * h;

    return points;
}

std::vector<double> getLagrangeBasisCoeffs(
    const std::vector<double>& nodes,
    double funcValue,
    int i
) {
    int n = nodes.size();
    std::vector<double> coeffs(n, 0.0);
    coeffs[0] = 1.0;
    double denom = 1.0;
    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        denom *= (nodes[i] - nodes[j]);
        for (int k = n - 2; k >= 0; --k) {
            coeffs[k + 1] += coeffs[k];
            coeffs[k] *= -nodes[j];
        }
    }
    double multiply = funcValue / denom;
    for (double& c : coeffs) {
        c *= multiply;
    }
    return coeffs;
}

std::vector<double>
getFunctionValues(
    std::vector<double>& points,
    std::function<double(double)> func
) {
    std::vector<double> values;
    values.reserve(points.size());
    for (auto point : points) {
        values.push_back(func(point));
    }
    return values;
}

double interpolationFunction(double x) {
    double result = interpolationValues[0];
    for (size_t i = 1; i < interpolationValues.size(); ++i) {
        double xPow = 1;
        for (int j = 0; j < i; ++j) {
            xPow *= x;
        }
        result += xPow * interpolationValues[i];
    }
    return result;
}


double interpolationErrorFunction(double x) {
    return std::abs(interpolationFunction(x) - function(x));
}


std::vector<double> getInterpolationValues(
    double xMin,
    double xMax,
    int valuesCount,
    std::function<double(double)> func
) {
    std::vector<double> values(valuesCount, 0.0);

    // std::vector<double> nodes =
    //     getEquidistantPoints(xMin, xMax, valuesCount);
    std::vector<double> nodes =
        getChebyshevNodes(valuesCount);

    std::vector<double> funcValues =
        getFunctionValues(nodes, func);

    for (int i = 0; i < valuesCount; ++i) {
        std::vector<double> coeffs =
            getLagrangeBasisCoeffs(
                nodes,
                funcValues[i],
                i
            );

        for (int j = 0; j < valuesCount; ++j) {
            values[j] += coeffs[j];
        }
    }

    return values;
}



int main() {
    double xMin, xMax;
    int valuesCount;
    std::cin >> xMin >> xMax >> valuesCount;
    interpolationValues = getInterpolationValues(xMin, xMax, valuesCount, function);
    std::vector<std::function<double(double)>> funcs = {function, interpolationFunction};
    std::vector<std::function<double(double)>> errFunc = {interpolationErrorFunction};
    std::vector<std::string> funcsName = {"Функция", "Интерполяция"};
    std::vector<std::string> errFuncName = {"|f(x) - P_n(x)|"};
    buildGraphs(xMin, xMax, funcs, funcsName, "testCheb.svg");
    buildGraphs(xMin, xMax, errFunc, errFuncName, "errCheb.svg");
}