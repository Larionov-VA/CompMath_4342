#include "./functions.hpp"
#include "./Newton-Simpson.hpp"
#include "./rectangle.hpp"
#include "./trapezoid.hpp"


void printListOfFunctions() {
    std::cout << "Выберете функцию:\n";
    std::cout << "1. Линейная\n2. Квадратичная\n3. Степенная\n";
    std::cout << "4. Показательная\n5. Логарифмическая\n";
    std::cout << "6. Синус\n7. Косинус\n8. Тангенс\n9. Котангенс\n";
}


void printListOfMethods() {
    std::cout << "Выберете метод вычисления интеграла:\n";
    std::cout << "1. Метод левых прямоугольников\n";
    std::cout << "2. Метод правых прямоугольников\n";
    std::cout << "3. Метод центральных прямоугольников\n";
    std::cout << "4. Метод трапеций\n";
    std::cout << "5. Метод Ньютона-Симпсона\n";
    std::cout << "0. Сравнение методов\n";
}


std::function<ld(ld)> getFunction(int functionChoice) {
    std::function<ld(ld)> func;
    switch (functionChoice) {
    case LINEAR:
        func = linear;
        std::cout << "Вы выбрали линейную функцию, она задается формулой:\n";
        std::cout << "\tf(x) = a*x + b\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case QUAD:
        func = quadratic;
        std::cout << "Вы выбрали квадратичную функцию, она задается формулой:\n";
        std::cout << "\tf(x) = a*x^2 + b*x + c\nВведите коэффициенты a, b и c.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b >> funcCoeff.c;
        break;
    case POW:
        func = power;
        std::cout << "Вы выбрали степенную функцию, она задается формулой:\n";
        std::cout << "\tf(x) = a*x^b\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case EXP:
        func = exponential;
        std::cout << "Вы выбрали показательную функцию, она задается формулой:\n";
        std::cout << "\ta*b^x\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case LOG:
        func = lg;
        std::cout << "Вы выбрали логарифмическую функцию, она задается формулой:\n";
        std::cout << "\ta + log2(x + b)\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case SIN:
        func = sn;
        std::cout << "Вы выбрали функцию синуса, она задается формулой:\n";
        std::cout << "\ta + sin(x + b)\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case COS:
        func = cs;
        std::cout << "Вы выбрали функцию косинуса, она задается формулой:\n";
        std::cout << "\ta + cos(x + b)\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case TG:
        func = tg;
        std::cout << "Вы выбрали функцию тангенса, она задается формулой:\n";
        std::cout << "\ta + tg(x + b)\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    case CTG:
        func = cg;
        std::cout << "Вы выбрали функцию котангенса, она задается формулой:\n";
        std::cout << "\ta + ctg(x + b)\nВведите коэффициенты a и b.\n";
        std::cin >> funcCoeff.a >> funcCoeff.b;
        break;
    default:
        func = linear;
        std::cout << "Выбрана функия по умолчанию:\n";
        std::cout << "\tf(x) = x\n";
        funcCoeff.a = 1;
        funcCoeff.b = 0;
        break;
    }
    return func;
}


std::function<ld(std::function<ld(ld)>, ld, ld, int)> getMethod(int methodChoice) {
    std::function<ld(std::function<ld(ld)>, ld, ld, int)> method;
    switch (methodChoice) {
    case RECTANGLE + rectangleType::LEFT:
        method = rectangle;
        rectMethodType = rectangleType::LEFT;
        std::cout << "Вы выбрали метод левых прямоугольников\n";
        break;
    case RECTANGLE + rectangleType::RIGHT:
        method = rectangle;
        rectMethodType = rectangleType::RIGHT;
        std::cout << "Вы выбрали метод правых прямоугольников\n";
        break;
    case RECTANGLE + rectangleType::CENTER:
        method = rectangle;
        rectMethodType = rectangleType::CENTER;
        std::cout << "Вы выбрали метод центральных прямоугольников\n";
        break;
    case TRAPEZOID:
        method = trapezoid;
        std::cout << "Вы выбрали метод трапеций\n";
        break;
    case SIMPSON:
        method = NewtonSimpson;
        std::cout << "Вы выбрали метод Ньютона-Симпсона\n";
        break;
    default:
        method = NewtonSimpson;
        std::cout << "метод Ньютона-Симпсона выбран по умолчанию\n";
        break;
    }
    return method;
}


int main() {
    printListOfFunctions();
    int functionChoice;
    std::cin >> functionChoice;
    auto function = getFunction(functionChoice);
    printListOfMethods();
    int methodChoice;
    std::cin >> methodChoice;
    if (methodChoice) {
        auto method = getMethod(methodChoice);
        int n, a, b;
        std::cout << "Введите пределы интегрирования.\n";
        std::cin >> a >> b;
        std::cout << "Введите количество отрезков.\n";
        std::cin >> n;
        std::cout << method(function, a, b, n) << '\n';
    }
    else {
        int n, a, b;
        std::cout << "Введите пределы интегрирования.\n";
        std::cin >> a >> b;
        std::cout << "Введите количество отрезков.\n";
        std::cin >> n;
        rectMethodType = rectangleType::LEFT;
        std::cout << "1. ";
        std::cout << rectangle(function, a, b, n) << '\n';
        rectMethodType = rectangleType::RIGHT;
        std::cout << "2. ";
        std::cout << rectangle(function, a, b, n) << '\n';
        rectMethodType = rectangleType::CENTER;
        std::cout << "3. ";
        std::cout << rectangle(function, a, b, n) << '\n';
        std::cout << "4. ";
        std::cout << trapezoid(function, a, b, n) << '\n';
        std::cout << "5. ";
        std::cout << NewtonSimpson(function, a, b, n) << '\n';
    }
    return EXIT_SUCCESS;
}