- 可调整算法：

  - 平流过程中采样不在网格内，目前采用最近点，后续可替换为基于水平集的外推等方法：

    ```c++
    inline double GetAt(const ScalarArray2D& d, int w, int h, int x, int y) {
      x = std::max(0, std::min(x, w - 1));
      y = std::max(0, std::min(y, h - 1));
      return d(x, y);
    }

    double SampleNearest(const ScalarArray2D& data, int width, int height, double x, double y) {
      int xi = static_cast<int>(std::round(Clamp(x, 0.0, static_cast<double>(width - 1))));
      int yi = static_cast<int>(std::round(Clamp(y, 0.0, static_cast<double>(height - 1))));
      return GetAt(data, width, height, xi, yi);
    }

    ......
    ```

- 需要注意的参数：
