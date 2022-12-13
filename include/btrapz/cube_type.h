struct Cube;
struct Cube
{
	int beg_t, end_t;
	double t;
	double t_dif;
	double beg_l, end_l;
	double upp_skew, upp_bias, down_skew, down_bias;
	double l_upp_skew, l_upp_bias, l_down_skew, l_down_bias;
	bool merge, split;
	int count;
	Cube()
	{		
		t = 0.0;
		beg_t = 0; end_t = 0;
		beg_l=0.0, end_l=0.0;
		upp_skew = 0.0; upp_bias = 1000.0; down_skew = 0.0; down_bias = 0.0;
		l_upp_skew = 0.0; l_upp_bias = 1000.0; l_down_skew = 0.0; l_down_bias = 0.0;
		merge = 0; split = 0;
		count = 0;
	}

	~Cube(){}
};
