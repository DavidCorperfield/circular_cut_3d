#ifndef _ROTATE_3D_H
#define _ROTATE_3D_H

#include "tpoint.h"
#include <vector>

using namespace std;

namespace ROTATE_3D {

	template<typename TYPE>
	class RotateMarix {
	public:
		vector<vector<TYPE>> R;  //Ðý×ª¾ØÕóR 3*3 °´ÐÐ´æ´¢
		vector<vector<TYPE>> RT;

		RotateMarix() { 
			R.resize(3);
			for (int i = 0; i < 3; i++)
			{
				R[i].resize(3);
			}

			RT.resize(3);
			for (int i = 0; i < 3; i++)
			{
				RT[i].resize(3);
			}
		}

		void GeneratR(TYPE cosPhi, TYPE sinPhi, TYPE cosOmiga, TYPE sinOmiga, TYPE cosKapa, TYPE sinKapa) {
			R[0][0] = cosPhi * cosKapa - sinPhi * sinOmiga * sinKapa;
			R[0][1] = - cosPhi * cosKapa - sinPhi * sinOmiga * sinKapa;
			R[0][2] = - sinPhi * cosOmiga;
			R[1][0] = cosOmiga * sinKapa;
			R[1][1] = cosOmiga * cosKapa;
			R[1][2] = - sinOmiga;
			R[2][0] = sinPhi * cosKapa + cosPhi * sinOmiga * sinKapa;
			R[2][1] = - sinPhi * sinKapa + cosPhi * sinOmiga * cosKapa;
			R[2][2] = cosPhi * cosOmiga;

			RT[0][0] = R[0][0];
			RT[0][1] = R[1][0];
			RT[0][2] = R[2][0];
			RT[1][0] = R[0][1];
			RT[1][1] = R[1][1];
			RT[1][2] = R[2][1];
			RT[2][0] = R[0][2];
			RT[2][1] = R[1][2];
			RT[2][2] = R[2][2];
		}

		template<typename TY_>
		void RotateR(TY_ ix, TY_ iy, TY_ iz, TY_ &ox, TY_ &oy, TY_ &oz)
		{
			ox = TY_(R[0][0] * ix + R[0][1] * iy + R[0][2] * iz);
			oy = TY_(R[1][0] * ix + R[1][1] * iy + R[1][2] * iz);
			oz = TY_(R[2][0] * ix + R[2][1] * iy + R[2][2] * iz);
		}

		template<typename TY_>
		void RotateRT(TY_ ix, TY_ iy, TY_ iz, TY_ &ox, TY_ &oy, TY_ &oz)
		{
			ox = TY_(RT[0][0] * ix + RT[0][1] * iy + RT[0][2] * iz);
			oy = TY_(RT[1][0] * ix + RT[1][1] * iy + RT[1][2] * iz);
			oz = TY_(RT[2][0] * ix + RT[2][1] * iy + RT[2][2] * iz);
		}
	};

} // namespace ROTATE_3D

#endif // _ROTATE_3D_H