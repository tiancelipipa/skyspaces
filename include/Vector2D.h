#pragma once

#include <cmath>

namespace skyspaces {

template <typename T>
struct Vector2 {
    T x{};
    T y{};

    Vector2() noexcept = default;
    Vector2(T x_value, T y_value) noexcept
        : x(x_value), y(y_value) {}

    Vector2& operator+=(const Vector2& other) noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2& operator-=(const Vector2& other) noexcept {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2& operator*=(T scalar) noexcept {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vector2& operator/=(T scalar) noexcept {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    friend Vector2 operator+(Vector2 lhs, const Vector2& rhs) noexcept {
        lhs += rhs;
        return lhs;
    }

    friend Vector2 operator-(Vector2 lhs, const Vector2& rhs) noexcept {
        lhs -= rhs;
        return lhs;
    }

    friend Vector2 operator*(Vector2 lhs, T scalar) noexcept {
        lhs *= scalar;
        return lhs;
    }

    friend Vector2 operator*(T scalar, Vector2 rhs) noexcept {
        rhs *= scalar;
        return rhs;
    }

    friend Vector2 operator/(Vector2 lhs, T scalar) noexcept {
        lhs /= scalar;
        return lhs;
    }

    friend Vector2 operator-(Vector2 vec) noexcept {
        vec.x = -vec.x;
        vec.y = -vec.y;
        return vec;
    }

    friend T Dot(const Vector2& a, const Vector2& b) noexcept {
        return a.x * b.x + a.y * b.y;
    }

    friend T Cross(const Vector2& a, const Vector2& b) noexcept {
        return a.x * b.y - a.y * b.x;
    }

    friend Vector2 Lerp(const Vector2& a, const Vector2& b, T t) noexcept {
        return a * (static_cast<T>(1) - t) + b * t;
    }

    T Length() const noexcept {
        return std::sqrt(x * x + y * y);
    }

    T LengthSquared() const noexcept {
        return x * x + y * y;
    }

    Vector2 Normalized() const noexcept {
        const T len = Length();
        return len > static_cast<T>(0)
            ? Vector2(x / len, y / len)
            : Vector2(static_cast<T>(0), static_cast<T>(0));
    }
};

using Vector2D = Vector2<double>;
using Vector2F = Vector2<float>;
using Vector2I = Vector2<int>;

} // namespace skyspaces
