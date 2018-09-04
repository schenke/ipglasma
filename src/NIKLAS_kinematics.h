#ifndef kinematics_h
#define kinematics_h



namespace KINEMATICS
{
    
    void SGL_to_CTR(DOUBLE x1, DOUBLE y1, DOUBLE x2, DOUBLE y2, DOUBLE &bx, DOUBLE &by, DOUBLE &rx, DOUBLE &ry)
    {
        bx = 0.5*(x1+x2);
        by = 0.5*(y1+y2);
        rx = x1-x2;
        ry = y1-y2;
    }
    
    
    void CTR_to_SGL(DOUBLE bx, DOUBLE by, DOUBLE rx, DOUBLE ry, DOUBLE &x1, DOUBLE &y1, DOUBLE &x2, DOUBLE &y2)
    {
        x1 = bx+0.5*rx;
        y1 = by+0.5*ry;
        x2 = bx-0.5*rx;
        y2 = by-0.5*ry;
    }

    
    
    
}




#endif
