#ifndef COLORMAPS_H
#define COLORMAPS_H
#include <cmath>
#include <tuple>


using namespace std;

class ColorMaps {
public:
    ColorMaps();
    virtual ~ColorMaps();
    static tuple<unsigned char,unsigned char,unsigned char> greyScale(const double);
    static tuple<unsigned char,unsigned char,unsigned char> viridis(const double);
    static tuple<unsigned char,unsigned char,unsigned char> plasma(const double);
    static tuple<unsigned char,unsigned char,unsigned char> inferno(const double);
    static tuple<unsigned char,unsigned char,unsigned char> magma(const double);
    static tuple< unsigned char,unsigned char,unsigned char> cividis(const double);
    static tuple<unsigned char,unsigned char,unsigned char> jet(const double);

protected:

private:
};

#endif // COLORMAPS_H
