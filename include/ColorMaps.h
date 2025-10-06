#ifndef COLORMAPS_H
#define COLORMAPS_H
#include <cmath>
#include <tuple>


using namespace std;

class ColorMaps
{
    public:
        ColorMaps();
        virtual ~ColorMaps();
        static tuple<unsigned char,unsigned char,unsigned char> greyScale(double);
        static tuple<unsigned char,unsigned char,unsigned char> viridis(double);
        static tuple<unsigned char,unsigned char,unsigned char> plasma(double);
        static tuple<unsigned char,unsigned char,unsigned char> inferno(double);
        static tuple<unsigned char,unsigned char,unsigned char> magma(double);
        static tuple< unsigned char,unsigned char,unsigned char> cividis(double);
        static tuple<unsigned char,unsigned char,unsigned char> jet(double);

    protected:

    private:
};

#endif // COLORMAPS_H
