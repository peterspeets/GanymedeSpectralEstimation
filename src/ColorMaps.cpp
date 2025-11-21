#include "ColorMaps.h"

ColorMaps::ColorMaps() {
    //ctor
}

tuple<unsigned char,unsigned char,unsigned char> ColorMaps::greyScale(const double pixelValue) {
    /*
     pixelValue is assumed to be between 0 and 1.
    */

    tuple<unsigned char,unsigned char,unsigned char> RGBValue = {static_cast<unsigned char>(255*pixelValue),static_cast<unsigned char>(255*pixelValue),static_cast<unsigned char>(255*pixelValue)};
    return RGBValue;
}

tuple<unsigned char,unsigned char,unsigned char> ColorMaps::viridis(const double pixelValue) {
    /*
    pixelValue is assumed to be between 0 and 1.
    */
    tuple<unsigned char,unsigned char,unsigned char> RGBValue;
    std::get<0>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.12477873223415754 - 0.09412643961695547 * pixelValue + 0.666877588575673 *pow(pixelValue,2) +
                                      1.159776225703296 *pow(pixelValue,3) + 0.23302421874297347*pow(pixelValue,4) -
                                      1.221574733565396 *pow(pixelValue,5) - 0.7394033860569665*pow(pixelValue,6) +
                                      0.5284442206334363 *pow(pixelValue,7) + 0.34897751079004824 *pow(pixelValue,8)));
    std::get<1>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.5650939419573958 + 0.46859984684845307 *pixelValue + 0.00433041277149867 *pow(pixelValue,2)+
                                      0.011504695781715193 *pow(pixelValue,3) - 0.23875005303171923*pow(pixelValue,4)-
                                      0.08044179641583385 *pow(pixelValue,5) + 0.17877730730557587 *pow(pixelValue,6) +
                                      0.05302885643416455 *pow(pixelValue,7) - 0.055987236950921386 *pow(pixelValue,8)));
    std::get<2>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.5510825822502813 - 0.09824202577332701 * pixelValue - 0.33962770156632033 *pow(pixelValue,2) -
                                      0.24186479476256878 *pow(pixelValue,3) + 0.1323435164578478 *pow(pixelValue,4) -
                                      0.03782652641530488 *pow(pixelValue,5) - 0.7000559918501962 *pow(pixelValue,6) +
                                      0.2876102051318723 *pow(pixelValue,7) + 0.5976866599222113*pow(pixelValue,8)));
    return RGBValue;
}


tuple<unsigned char,unsigned char,unsigned char> ColorMaps::plasma(const double pixelValue) {
    /*
    pixelValue is assumed to be between 0 and 1.
    */
    tuple<unsigned char,unsigned char,unsigned char> RGBValue;
    std::get<0>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.7964528741101662 + 0.4675081381727712 *pixelValue - 0.23019996550150998 *pow(pixelValue,2) +
                                      0.1286054864893778 *pow(pixelValue,3) - 0.15080463534143457 *pow(pixelValue,4) -
                                      0.35458431117789596 *pow(pixelValue,5) + 0.23341273319976882 *pow(pixelValue,6) +
                                      0.2014357906289061 *pow(pixelValue,7) - 0.15523349503440673 *pow(pixelValue,8)));
    std::get<1>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.2802581760154954 + 0.5844303201561526 *pixelValue - 0.10512461141927082 *pow(pixelValue,2) -
                                      0.013843939635487877 *pow(pixelValue,3) + 0.9474194884309068 *pow(pixelValue,4) -
                                      0.36843132864367617 *pow(pixelValue,5) - 0.9150299918959578 *pow(pixelValue,6) +
                                      0.27975914813246 *pow(pixelValue,7) + 0.2942978448041723 *pow(pixelValue,8)));
    std::get<2>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.47053105290139957 - 0.4636327077378959 *pixelValue + 0.13447654325019023 *pow(pixelValue,2) +
                                      0.19470573329359556 *pow(pixelValue,3) - 1.1065188378066797 *pow(pixelValue,4) +
                                      0.17947402624298264 *pow(pixelValue,5) + 1.4671087319228278 *pow(pixelValue,6) -
                                      0.10003120037060716 *pow(pixelValue,7) - 0.6297623507353992 *pow(pixelValue,8)));
    return RGBValue;
}


tuple<unsigned char,unsigned char,unsigned char> ColorMaps::inferno(const double pixelValue) {
    /*
    pixelValue is assumed to be between 0 and 1.
    */
    tuple<unsigned char,unsigned char,unsigned char> RGBValue;
    std::get<0>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.7330043366913852 + 0.7259448283222094 *pixelValue - 0.2921934546431462 *pow(pixelValue,2) -
                                      0.36534603625606077 *pow(pixelValue,3) + 0.06473346453583473 *pow(pixelValue,4) +
                                      0.13085459384570364 *pow(pixelValue,5) - 0.4745998536816382 *pow(pixelValue,6) -
                                      0.0027352197616549987 *pow(pixelValue,7) + 0.4670039951077641 *pow(pixelValue,8)));
    std::get<1>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.21610651636152473 + 0.41277725954148875 *pixelValue + 0.4059883587948613 *pow(pixelValue,2) +
                                      0.46458235861242964 *pow(pixelValue,3) - 0.23078259051210623 *pow(pixelValue,4) -
                                      0.7473217647021523 *pow(pixelValue,5) + 0.43857925555650695 *pow(pixelValue,6) +
                                      0.38106487380062215 *pow(pixelValue,7) - 0.33844568536951364 *pow(pixelValue,8)));
    std::get<2>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.33624224097455696 - 0.4245987913018207 *pixelValue - 0.611384188290079 *pow(pixelValue,2) -
                                      0.43647385641613906 *pow(pixelValue,3) + 1.0539870151090227 *pow(pixelValue,4) +
                                      2.843486831505292 *pow(pixelValue,5) - 0.5460585005651183 *pow(pixelValue,6) -
                                      1.6743800445300523 *pow(pixelValue,7) + 0.09522055226736005 *pow(pixelValue,8)));
    return RGBValue;
}

tuple<unsigned char,unsigned char,unsigned char> ColorMaps::magma(const double pixelValue) {
    /*
    pixelValue is assumed to be between 0 and 1.
    */
    tuple<unsigned char,unsigned char,unsigned char> RGBValue;
    std::get<0>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.7142969649682965 + 0.8214552410735821 *pixelValue - 0.10304747572792063 *pow(pixelValue,2) -
                                      0.7604619555349517 *pow(pixelValue,3) - 0.9058561146957473 *pow(pixelValue,4) +
                                      0.8497534964270057 *pow(pixelValue,5) + 1.3461142497817657 *pow(pixelValue,6) -
                                      0.4225315305120873 *pow(pixelValue,7) - 0.5611147694432717 *pow(pixelValue,8)));
    std::get<1>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.2141216109987368 + 0.3113458772800889 *pixelValue + 0.24796291244569338 *pow(pixelValue,2) +
                                      0.8300668152855053 *pow(pixelValue,3) + 0.6911637602542952 *pow(pixelValue,4) -
                                      1.3644841451209861 *pow(pixelValue,5) - 1.203114782467431 *pow(pixelValue,6) +
                                      0.7368272097958986 *pow(pixelValue,7) + 0.5455731566039694 *pow(pixelValue,8)));
    std::get<2>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.47319250768070903 - 0.28857966234699983 *pixelValue - 0.4357116634863791 *pow(pixelValue,2) +
                                      0.6458485053938033 *pow(pixelValue,3) + 1.5881642745860713 *pow(pixelValue,4) +
                                      0.6419689985437824 *pow(pixelValue,5) - 2.544673110967218 *pow(pixelValue,6) -
                                      0.6406906128472649 *pow(pixelValue,7) + 1.3133148115774373 *pow(pixelValue,8)));
    return RGBValue;
}



tuple<unsigned char,unsigned char,unsigned char> ColorMaps::cividis(const double pixelValue) {
    /*
    pixelValue is assumed to be between 0 and 1.
    */
    tuple<unsigned char,unsigned char,unsigned char> RGBValue;
    std::get<0>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.48780670161886697 + 0.4690174751636906 *pixelValue + 0.08115097657044208 *pow(pixelValue,2) -
                                      0.1743169832292792 *pow(pixelValue,3) - 0.044122319603577306 *pow(pixelValue,4) +
                                      0.8162169383691712 *pow(pixelValue,5) - 0.4909717782172515 *pow(pixelValue,6) -
                                      0.6191803659213162 *pow(pixelValue,7) + 0.4814174325337144 *pow(pixelValue,8)));
    std::get<1>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.48381046811949757 + 0.3663352027961549 *pixelValue + 0.03656326050064862 *pow(pixelValue,2) +
                                      0.029339858700067913 *pow(pixelValue,3) - 0.004735610494240279 *pow(pixelValue,4) -
                                      0.020691456158095246 *pow(pixelValue,5) + 0.0005710088286790926 *pow(pixelValue,6) +
                                      0.011181286587903039 *pow(pixelValue,7) + 0.005498983392409244 *pow(pixelValue,8)));
    std::get<2>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.46770210642764803 + 0.07170924197990428 *pixelValue - 0.3204713077164482 *pow(pixelValue,2) -
                                      0.2580258658496583 *pow(pixelValue,3) + 0.7919931807546625 *pow(pixelValue,4) -
                                      0.10103077990190681 *pow(pixelValue,5) - 0.9841091183208727 *pow(pixelValue,6) +
                                      0.24235698598317848 *pow(pixelValue,7) + 0.2864874198456464 *pow(pixelValue,8)));
    return RGBValue;
}

tuple<unsigned char,unsigned char,unsigned char> ColorMaps::jet(const double pixelValue) {
    /*
    pixelValue is assumed to be between 0 and 1.
    */
    tuple<unsigned char,unsigned char,unsigned char> RGBValue;
    std::get<0>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.4854335860904151 + 1.5821094808402967 *pixelValue - 0.2963546645545918 *pow(pixelValue,2) +
                                      4.301902675716216 *pow(pixelValue,3) + 7.5062472041604975 *pow(pixelValue,4) -
                                      78.34623210466732 *pow(pixelValue,5) - 51.97200528298831 *pow(pixelValue,6) + 347.8438515993157 *pow(pixelValue,7) +
                                      164.04461629972764 *pow(pixelValue,8) - 747.4618178807872 *pow(pixelValue,9) -
                                      260.6521907328631 *pow(pixelValue,10) + 855.4362942742106 *pow(pixelValue,11) +
                                      202.9896067838633 *pow(pixelValue,12) - 502.9990340866969 *pow(pixelValue,13) -
                                      64.19399446386244 *pow(pixelValue,14) + 119.89663855781767 *pow(pixelValue,15) +
                                      2.3295470537353635 *pow(pixelValue,16)
                                     ));
    std::get<1>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.988186824872677 - 0.09906781565329965 *pixelValue + 2.2078947980580437 *pow(pixelValue,2) +
                                      4.3166938103288945 *pow(pixelValue,3) - 52.11912133428867 *pow(pixelValue,4) -
                                      35.57292450123659 *pow(pixelValue,5) + 274.40444635725805 *pow(pixelValue,6) +
                                      141.5922135365645 *pow(pixelValue,7) - 757.1797486865643 *pow(pixelValue,8) - 297.5700897044173 *pow(pixelValue,9) +
                                      1180.641252807954 *pow(pixelValue,10) + 333.86376201915743 *pow(pixelValue,11) -
                                      1029.3772929347822 *pow(pixelValue,12) - 188.3313263827873 *pow(pixelValue,13) +
                                      462.37864594585943 *pow(pixelValue,14) + 41.792786085942396 *pow(pixelValue,15) -
                                      81.93712375369341 *pow(pixelValue,16)));
    std::get<2>(RGBValue) =  static_cast<unsigned char>(
                                 255*(0.4854335860903685 - 1.5821094808403908 *pixelValue - 0.2963546645527013 *pow(pixelValue,2) -
                                      4.301902675710786 *pow(pixelValue,3) + 7.506247204145384 *pow(pixelValue,4) + 78.34623210459425 *pow(pixelValue,5) -
                                      51.97200528297464 *pow(pixelValue,6) - 347.84385159892884 *pow(pixelValue,7) +
                                      164.0446162999536 *pow(pixelValue,8) + 747.4618178797654 *pow(pixelValue,9) -
                                      260.65219073376704 *pow(pixelValue,10) - 855.4362942727736 *pow(pixelValue,11) +
                                      202.98960678536235 *pow(pixelValue,12) + 502.99903408566803 *pow(pixelValue,13) -
                                      64.1939944650587 *pow(pixelValue,14) - 119.89663855752293 *pow(pixelValue,15) + 2.32954705411107 *pow(pixelValue,16)));
    return RGBValue;
}

ColorMaps::~ColorMaps() {
    //dtor
}
