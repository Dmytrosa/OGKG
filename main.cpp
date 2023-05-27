#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <SFML/Graphics.hpp>

struct Point {
    double x, y;
};

struct Circle {
    Point center;
    double radius;
};

const double INF = std::numeric_limits<double>::infinity();

double crossProduct(const Point& p, const Point& q, const Point& r) {
    return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
}

bool comparePoints(const Point& p1, const Point& p2) {
    if (p1.x != p2.x)
        return p1.x < p2.x;
    return p1.y < p2.y;
}

std::vector<Point> findConvexHull(std::vector<Point>& points) {
    std::vector<Point> hull;

    if (points.size() <= 3) {
        return points;
    }

    std::sort(points.begin(), points.end(), comparePoints);

    std::vector<Point> upperHull;
    upperHull.push_back(points[0]);
    upperHull.push_back(points[1]);

    for (size_t i = 2; i < points.size(); ++i) {
        upperHull.push_back(points[i]);

        while (upperHull.size() >= 3 &&
            crossProduct(upperHull[upperHull.size() - 3], upperHull[upperHull.size() - 2], upperHull[upperHull.size() - 1]) <= 0) {
            upperHull.erase(upperHull.end() - 2);
        }
    }

    std::vector<Point> lowerHull;
    lowerHull.push_back(points.back());
    lowerHull.push_back(points[points.size() - 2]);

    for (int i = points.size() - 3; i >= 0; --i) {
        lowerHull.push_back(points[i]);

        while (lowerHull.size() >= 3 &&
            crossProduct(lowerHull[lowerHull.size() - 3], lowerHull[lowerHull.size() - 2], lowerHull[lowerHull.size() - 1]) <= 0) {
            lowerHull.erase(lowerHull.end() - 2);
        }
    }

    upperHull.pop_back();
    lowerHull.pop_back();

    hull.reserve(upperHull.size() + lowerHull.size());
    hull.insert(hull.end(), upperHull.begin(), upperHull.end());
    hull.insert(hull.end(), lowerHull.begin(), lowerHull.end());

    return hull;
}

Circle findLargestEmptyCircle(const std::vector<Point>& points, const std::vector<Point>& boundary) {
    Circle largestCircle;

    if (boundary.size() < 3) {
        double totalX = 0.0, totalY = 0.0;
        for (const auto& point : boundary) {
            totalX += point.x;
            totalY += point.y;
        }
        largestCircle.center = { totalX / boundary.size(), totalY / boundary.size() };
        largestCircle.radius = 0;
        return largestCircle;
    }

    const Point& p1 = boundary[0];
    const Point& p2 = boundary[1];
    const Point& p3 = boundary[2];

    double a = p2.x - p1.x;
    double b = p2.y - p1.y;
    double c = p3.x - p1.x;
    double d = p3.y - p1.y;
    double e = a * (p1.x + p2.x) + b * (p1.y + p2.y);
    double f = c * (p1.x + p3.x) + d * (p1.y + p3.y);
    double g = 2.0 * (a * (p3.y - p2.y) - b * (p3.x - p2.x));

    if (g == 0) {
        largestCircle.center = { 0.0, 0.0 };
        largestCircle.radius = -INF;
        return largestCircle;
    }

    double centerX = (d * e - b * f) / g;
    double centerY = (a * f - c * e) / g;

    largestCircle.center = { centerX, centerY };
    largestCircle.radius = std::hypot(p1.x - centerX, p1.y - centerY);


    for (const Point& point : points) {
        if (std::hypot(point.x - largestCircle.center.x, point.y - largestCircle.center.y) > largestCircle.radius) {
            std::vector<Point> newBoundary = boundary;
            newBoundary.push_back(point);

            Circle newCircle = findLargestEmptyCircle(points, newBoundary);

            if (newCircle.radius > largestCircle.radius) {
                largestCircle = newCircle;
            }
        }
    }

    double min_distance = std::numeric_limits<double>::max();
    for (size_t i = 0; i < boundary.size(); ++i) {
        const Point& p1 = boundary[i];
        const Point& p2 = boundary[(i + 1) % boundary.size()];

        double numerator = std::abs((p2.y - p1.y) * centerX - (p2.x - p1.x) * centerY + p2.x * p1.y - p2.y * p1.x);
        double denominator = std::hypot(p2.y - p1.y, p2.x - p1.x);

        double distance = numerator / denominator;
        if (distance < min_distance) {
            min_distance = distance;
        }
    }
    largestCircle.radius = min_distance;


    return largestCircle;
}

int main() {
    sf::RenderWindow window(sf::VideoMode(800, 800), "Largest Empty Circle");

    std::vector<Point> points = {
        {200.0, 200.0},
        {200.0, 0.0},
        {0.0, 0.0},
        {0.0, 200.0},
        {100.0, 300.0},
        {300.0, 300.0},
        {100.0, 100.0},
        {300.0, 100.0},
        {200.0, 100.0}
    };

    std::vector<Point> boundary = findConvexHull(points);
    Circle largestEmptyCircle = findLargestEmptyCircle(points, boundary);

    // Створюємо точки
    std::vector<sf::CircleShape> sfPoints;
    for (const Point& point : points) {
        sf::CircleShape sfPoint(5);
        sfPoint.setPosition(point.x, point.y);
        sfPoint.setFillColor(sf::Color::Red);
        sfPoints.push_back(sfPoint);
    }

    // Створюємо лінії
    std::vector<sf::RectangleShape> sfLines;
    for (size_t i = 0; i < boundary.size(); ++i) {
        const Point& p1 = boundary[i];
        const Point& p2 = boundary[(i + 1) % boundary.size()];

        double distance = std::hypot(p2.x - p1.x, p2.y - p1.y);
        sf::RectangleShape sfLine(sf::Vector2f(distance, 2));
        sfLine.setPosition(p1.x, p1.y);
        sfLine.setRotation(std::atan2(p2.y - p1.y, p2.x - p1.x) * 180 / 3.14159265);
        sfLine.setFillColor(sf::Color::Green);

        sfLines.push_back(sfLine);
    }

    // Створюємо коло
    sf::CircleShape sfCircle(largestEmptyCircle.radius);
    sfCircle.setOutlineThickness(3);
    sfCircle.setOutlineColor(sf::Color::White);
    sfCircle.setFillColor(sf::Color::Transparent);
    sfCircle.setPosition(largestEmptyCircle.center.x - largestEmptyCircle.radius, largestEmptyCircle.center.y - largestEmptyCircle.radius);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        window.clear();

        // Виводимо об'єкти
        for (const sf::CircleShape& sfPoint : sfPoints) {
            window.draw(sfPoint);
        }
        for (const sf::RectangleShape& sfLine : sfLines) {
            window.draw(sfLine);
        }
        window.draw(sfCircle);

        window.display();
    }

    return 0;
}