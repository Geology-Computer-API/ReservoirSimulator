//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#include "TRSRibFrac.h"

TRSRibFrac::TRSRibFrac(){
     faxis.Resize(3, 3);
     fdata.Resize(4, 3);
}

TRSRibFrac::TRSRibFrac(Matrix data){
     faxis.Resize(3, 3);
     Check_ConsistencyData(data);
}

void TRSRibFrac::SetPlane(Matrix plane){
     faxis.Resize(3, 3);
     Check_ConsistencyData(plane);
}

Matrix TRSRibFrac::GetPlane() const{
    return fdata;
}

void TRSRibFrac::SetTolerance(double tolerance){
    ftolerance = tolerance;
}

double TRSRibFrac::GetTolerance() const{
    return ftolerance;
}

/**
 * @brief Consistency of plane points
 * @param Matrix nx3, n is the number of points with the axis x, y , z
 * @return False if the points are not coplanar
 */

void TRSRibFrac::Check_ConsistencyData(Matrix data) {
  
// Checking vector consistency
    int cols = data.Cols();
    int rows = data.Rows();
    if(cols != 3){
        std::cout<<"Check the input data";
        DebugStop();
    }
    if(rows < 3 or rows > 4){
        std::cout<<"Check the input data";
        DebugStop();
    }
    
//Ax0 computation
    faxis(0,0)=data(1,0)-data(0,0);
    faxis(0,1)=data(1,1)-data(0,1);
    faxis(0,2)=data(1,2)-data(0,2);
    
//Ax1 without normalization
    faxis(1,0)=data(rows-1,0)-data(0,0);
    faxis(1,1)=data(rows-1,1)-data(0,1);
    faxis(1,2)=data(rows-1,2)-data(0,2);
    
//Ax1 normalization
    double normx = faxis(1,0) - (faxis(0,0))*((faxis(0,0)*faxis(1,0)) + (faxis(0,1)*faxis(1,1)) + (faxis(0,2)*faxis(1,2)));
    
    double normy = faxis(1,1)- faxis(0,1)*(faxis(0,0)*faxis(1,0) + faxis(0,1)*faxis(1,1) + faxis(0,2)*faxis(1,2));
    
    double normz = faxis(1,2)- faxis(0,2)*(faxis(0,0)*faxis(1,0) + faxis(0,1)*faxis(1,1) + faxis(0,2)*faxis(1,2));
    
    faxis(1,0)=normx;
    faxis(1,1)=normy;
    faxis(1,2)=normz;
    
//Ax2 computation
    faxis(2,0)=faxis(0,1)*faxis(1,2) - faxis(0,2)*faxis(1,1);
    faxis(2,1)=faxis(0,2)*faxis(1,0) - faxis(0,0)*faxis(1,2);
    faxis(2,2)=faxis(0,0)*faxis(1,1) - faxis(0,1)*faxis(1,0);

    if (rows ==4 && cols ==3){
        
//Coplanar verification
        double ver = faxis(2,0)*(data(2,0)-data(0,0))+faxis(2,1)*(data(2,1)-data(0,1))+faxis(2,2)*(data(2,2)-data(0,2));
        
        if(abs(ver) > ftolerance){
            std::cout<<"The points are not coplanar"<<"\n"<<std::endl;
            DebugStop();
        }
    }
    
//Mid points computation
    data.Resize(rows+rows, 3);
    for(int i=0; i<(rows); i++){
        data(rows+i,0)=(data(i,0)+data(i+1,0))/2;
        data(rows+i,1)=(data(i,1)+data(i+1,1))/2;
        data(rows+i,2)=(data(i,2)+data(i+1,2))/2;
    }
    data(rows+rows-1,0)=(data((rows)-1,0)+data(0,0))/2;
    data(rows+rows-1,1)=(data((rows)-1,1)+data(0,1))/2;
    data(rows+rows-1,2)=(data((rows)-1,2)+data(0,2))/2;
    
    data.Resize(rows+rows+1, 3);
    for(int i=0; i< rows; i++){
    
            data(rows+rows,0) += (1.0/rows)*data(i,0);
            data(rows+rows,1) += (1.0/rows)*data(i,1);
            data(rows+rows,2) += (1.0/rows)*data(i,2);
    }
    
    fdata = data;
}

/**
 * @brief Check if a point is above or below a plane
 * @param Point
 * @return True if the point is above the plane and false if is below the planealerta
 */

//Checking if the given point is above the plane
bool TRSRibFrac::Check_point_above(TPZVec<double> point) const{
    int pm = fdata.Rows();
    double point_distance = (point[0] - fdata(pm-1,0))*faxis(2,0) + (point[1] - fdata(pm-1,1))*faxis(2,1)+(point[2] - fdata(pm-1,2))*faxis(2,2);
    if (point_distance>0){
        return true;
    }
    else{
        return false;
    };
};

/**
 * @brief Check two points if they are on both sides of the plan
 * @param Point 1 and Point 2
 * @return True if both are on two sides of the plane, false otherwise
 */

bool TRSRibFrac::Check_rib(TPZVec<double> p1, TPZVec<double> p2) const {
    if(Check_point_above(p1)==true && Check_point_above(p2)==false){
        return true;
    }
    else if (Check_point_above(p1)==false && Check_point_above(p2)==true){
        return true;
    }
    else{
        return false;
    };
};

/**
 * @brief Check if the neighbour has a lower dimension
 * @param Geo element side
 * @return True if has a lower dimension
 */

bool TRSRibFrac::HasLowerDimensionNeighbour(TPZGeoElSide &gelside){
    int dimension = gelside.Dimension();
    if (gelside.Element()->Dimension() == dimension){
        return true;
    }

    TPZGeoElSide neighbour = gelside.Neighbour();
    
    while (neighbour != gelside){
        if (neighbour.Element()->Dimension()==dimension){
            return true;
            neighbour = neighbour.Neighbour();
         }
    return false;
    }
}

/**
 * @brief Creates the skeleton mesh
 * @param Dimension and number material ID
 * @return
 */

void TRSRibFrac::CreateSkeleton(int dimension, int matid){
    if(fmesh){
        int nel = fmesh->NElements();
        for(int iel=0; iel<nel; iel++){
            TPZGeoEl *gel = fmesh->Element(iel);
            int nsides = gel->NSides();
            for(int iside=0; iside<nsides; iside++){
                TPZGeoElSide gelside = gel->Neighbour(iside);
                if (gelside.Dimension()==dimension){
                    bool haskel = HasLowerDimensionNeighbour(gelside);
                    if(haskel==false){
                        TPZGeoElBC(gelside,matid);
                    }
                }
            }
        }
        
    }
    else{
        std::cout<<"";
    }
}
