#include "./includes.hpp"



std::vector<double> getChebyshevNodes(int nodesCount) {
    std::vector<double> nodes(nodesCount);
    for (int i = 0; i < nodesCount; ++i) {
        nodes[i] = std::cos((2.0L * (i + 1) - 1.0L) * M_PI / (2.0L * nodesCount));
    }
    return nodes;
}



void buildGraphs(
    const double lowerXbound,
    const double upperXbound,
    std::vector<std::function<double(double)>>& functions,
    std::vector<std::string>& functionsName,
    const std::string& filename
) {
    setenv("GNUTERM", "svg enhanced background rgb 'white'", 1);
    auto x = linspace(lowerXbound, upperXbound, std::abs(upperXbound - lowerXbound) * 100);
    figure(true);
    gcf()->color("white");
    std::vector<std::string> names(functions.size());
    for (size_t i = 0; i < functions.size(); ++i) {
        auto y = transform(x, functions[i]);
        auto p = plot(x, y)->line_width(1);
        hold(on);
    }
    hold(off);
    auto l = ::matplot::legend(functionsName);
    gca()->x_axis().zero_axis(true);
    gca()->y_axis().zero_axis(true);
    grid(on);
    save(filename);
}


double function(double x) {
    return std::pow(x, 3) - 3 * x - 2 * std::exp(-x);
}


double derivative1(double x) {
    return 2 * std::pow(M_El, -x) + 3 * std::pow(x, 2) - 3;
}


double derivative2(double x) {
    return - 2 * std::pow(M_El, -x) + 6 * x;
}


double derivative3(double x) {
    return 2 * std::pow(M_El, -x) + 6;
}


double derivative4(double x) {
    return - 2 * std::pow(M_El, -x);
}


double derivative5(double x) {
    return 2 * std::pow(M_El, -x);
}



