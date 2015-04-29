// FILL IN CODE (line below says "no" for all spheres, so replace it)
	Vect    w0, w;          // ray vectors
	float   a, b;           // params to calc ray-sphere intersect
	float t;
	Vect I;
	Intersection *inter;
	//float rayDeltaP = *(ray)->dir;
	Vect deltaP;
	VectSub(S->P, ray->orig, deltaP);
	float rayDeltaP = VectDotProd(ray->dir, deltaP);  //ray->dir * deltaP
	float tout = rayDeltaP; //ray->dir * deltaP;
	float tin = sqrt(pow((rayDeltaP), 2) - (pow(fabs(VectMag(deltaP)), 2) - pow(S->radius, 2)));
	// get triangle edge vectors and plane normal

	float t0 = tout + tin;
	float t1 = tout - tin;

	if(t0 > SMALL_NUM){			// ray intersects with sphere at t0
		if(t1 > SMALL_NUM){		// ray intersects with sphere at t1
			if(t0 < t1)			// if t0 < t1 then t0 is the first point
				t = t0;
			else				// else t1 is the first point
				t = t1;
		}
		else					// ray intersects with sphere at t0
			t = t0;
	}
	else if(t1 > SMALL_NUM)		// if t1 > SMALL_NUM ray intersects with sphere at t1 only
		t = t1;
	else						// else ray does not intersect with sphere at all or with positive t
		return NULL;



	// get intersect point of ray with triangle plane

	if (t < rayeps)                   // triangle is behind/too close to ray => no intersect
		return NULL;                 // for a segment, also test if (t > 1.0) => no intersect

	// intersect point of ray and plane

	VectAddS(t, ray->dir, ray->orig, I);

	inter = make_intersection();
	inter->t = t;
	VectCopy(inter->P, I);
	return inter;                      // I is in T
	//return NULL
