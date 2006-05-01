/** \file Galerkin.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface the class Galerkin.
 */
#ifndef GALERKIN_H_
#define GALERKIN_H_

#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "Stencil.h"
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

void computeGalerkin(
    NumericArray& resultL,
    PositionArray& resultJx,
    PositionArray& resultJy,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Restriction& restriction,
    const Stencil& stencil,
    const Prolongation& prolongation);

class Quadruple
{
private:
    Index v1_, v2_, v3_, v4_;
public:
    Quadruple(const Index v1,const Index v2,const Index v3,const Index v4)
        : v1_(v1),v2_(v2),v3_(v3),v4_(v4) {}

    Quadruple(const Quadruple& rhs)
        : v1_(rhs.v1_),v2_(rhs.v2_),v3_(rhs.v3_),v4_(rhs.v4_) {}
        
    Quadruple& operator=(const Quadruple& rhs)
    {
        v1_ = rhs.v1_;
        v2_ = rhs.v2_;
        v3_ = rhs.v3_;
        v4_ = rhs.v4_;
        return *this;
    }

    bool operator==(const Quadruple& rhs) const
    {
        if (this == &rhs)
            return true;
        return v1_==rhs.v1_ && v2_==rhs.v2_ && v3_== rhs.v3_ && v4_==rhs.v4_;
    }

    bool operator!=(const Quadruple& rhs) const
    {
        return !(*this==rhs);
    }
    bool operator<(const Quadruple& rhs) const
    {
        return 
               v1_<rhs.v1_ 
            || v1_==rhs.v1_ && v2_<rhs.v2_
            || v1_==rhs.v1_ && v2_==rhs.v2_ && v3_<rhs.v3_
            || v1_==rhs.v1_ && v2_==rhs.v2_ && v3_==rhs.v3_ && v4_<rhs.v4_;
    }
    bool operator>(const Quadruple& rhs) const
    {
        return (*this!=rhs) && !(*this<rhs);
    }
    bool operator>=(const Quadruple& rhs) const
    {
        return !(*this<rhs);
    }
    bool operator<=(const Quadruple& rhs) const
    {
        return !(*this>rhs);
    }
};

template<typename FineGridOperator>
class Galerkin : public Stencil
{
private:
    mutable std::vector<std::map<Quadruple,NumericArray > > data_;
    const FineGridOperator fineGridOperator_;
    Index size_;
    std::vector<std::vector<PositionArray > > jX_;
    std::vector<std::vector<PositionArray > > jY_;
    std::vector<const Prolongation*> prolongations_;
    std::vector<const Restriction*> restrictions_;
    Index currentDepth_;
    const Index nx_;
    const Index ny_;
    void update();
public:
    Galerkin(
        const FineGridOperator& fineGridOperator,
        const Index nx, const Index ny)
        : data_(9),
          fineGridOperator_(fineGridOperator),
          size_(fineGridOperator.size()),
          jX_(0),jY_(0),
          prolongations_(0),
          restrictions_(0),
          currentDepth_(prolongations_.size()),
          nx_(nx),
          ny_(ny){}
    virtual ~Galerkin() {}

    virtual Precision apply(
        const NumericArray& u,
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    virtual Precision getCenter(
      const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    virtual const NumericArray& getL(
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    inline const PositionArray& getJx(const Position p) const
    {
        if ( currentDepth_ == 0 )
            return fineGridOperator_.getJx( p );
        return jX_[currentDepth_][p];
    }

    inline const PositionArray& getJy(const Position p) const
    {
        if ( currentDepth_ == 0 )
            return fineGridOperator_.getJy( p );
        return jY_[currentDepth_][p];
    }

    void pushProlongation(const Prolongation& prolongation)
    {
        prolongations_.push_back(&prolongation);
        if ( prolongations_.size() == restrictions_.size() )
        {
            update();    
        }
    }

    void popProlongation()
    {
        prolongations_.pop_back();
        if ( prolongations_.size() == restrictions_.size() )
        {
            update();    
        }
    }

    void pushRestriction(const Restriction& restriction)
    {
        restrictions_.push_back(&restriction);
        if ( prolongations_.size() == restrictions_.size() )
        {
            update();    
        }
    }

    void popRestriction()
    {
        restrictions_.pop_back();
        if ( prolongations_.size() == restrictions_.size() )
        {
            update();    
        }
    }

    inline Index size() const
    {
        return size_;
    }

    inline bool isConstant() const
    {
        return fineGridOperator_.isConstant();
    }
};

template<typename FineGridOperator>
void Galerkin<FineGridOperator>::update()
{
    if ( jX_.size() >= prolongations_.size() )
    {
        currentDepth_ = prolongations_.size();
        return;
    }
    const Index divisor = static_cast<Index>(
                            std::pow( 2.0,
                                static_cast<Integer>( currentDepth_ + 1 ) ) );
    const Index minNxNy = 16;
    const Index nx = std::min( minNxNy, nx_/divisor );
    const Index ny = std::min( minNxNy, ny_/divisor );
    PositionArray jX;
    PositionArray jY;
    NumericArray operatorL;
    std::vector< PositionArray> jXvector;
    std::vector< PositionArray> jYvector;
    //C
    computeGalerkin(
        operatorL,jX,jY,
        C,nx/2,ny/2,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[C][Quadruple(nx/2,ny/2,nx,ny)]=operatorL;
    //W
    computeGalerkin(
        operatorL,jX,jY,
        W,1,ny/2,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[W][Quadruple(1,ny/2,nx,ny)]=operatorL;
    //N
    computeGalerkin(
        operatorL,jX,jY,
        N,nx/2,ny-1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[N][Quadruple(nx/2,ny-1,nx,ny)]=operatorL;
    //E
    computeGalerkin(
        operatorL,jX,jY,
        E,nx-1,ny/2,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[E][Quadruple(nx-1,ny/2,nx,ny)]=operatorL;
    //S
    computeGalerkin(
        operatorL,jX,jY,
        S,nx/2,1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[S][Quadruple(nx/2,1,nx,ny)]=operatorL;
    //NW
    computeGalerkin(
        operatorL,jX,jY,
        NW,1,ny-1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[NW][Quadruple(1,ny-1,nx,ny)]=operatorL;
    //NE
    computeGalerkin(
        operatorL,jX,jY,
        NE,nx-1,ny-1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[NE][Quadruple(nx-1,ny-1,nx,ny)]=operatorL;
    //SE
    computeGalerkin(
        operatorL,jX,jY,
        SE,nx-1,1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[SE][Quadruple(nx-1,1,nx,ny)]=operatorL;
    //SW
    computeGalerkin(
        operatorL,jX,jY,
        SW,1,1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    jXvector.push_back( jX );
    jYvector.push_back( jY );
    data_[SW][Quadruple(nx-1,1,nx,ny)]=operatorL;
    jX_.push_back(jXvector);
    jY_.push_back(jYvector);
    currentDepth_=prolongations_.size();
}

template<typename FineGridOperator>
Precision Galerkin<FineGridOperator>::apply(
  const NumericArray& u,
    const Position position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
  Precision result=0;
  NumericArray operatorL=getL(position,sx,sy,nx,ny);
  PositionArray jX=getJx(position);
  PositionArray jY=getJy(position);
  for (Index i=0; i<operatorL.size(); ++i)
      result+=operatorL[i]*u[(sy+jY[i])*(nx+1)+sx+jX[i]];
  return result;
}

template<typename FineGridOperator>
Precision Galerkin<FineGridOperator>::getCenter(
  const Position position,
  const Index sx,
  const Index sy,
  const Index nx,
  const Index ny) const
{
  return getL(position,sx,sy,nx,ny)[0];
}

template<typename FineGridOperator>
const NumericArray& Galerkin<FineGridOperator>::getL(
  const Position position,
  const Index sx,
  const Index sy,
  const Index nx,
  const Index ny) const
{
    static Index counter=0;
    ++counter;
    std::cout<<counter<<std::endl;
    if ( nx_ == nx && ny_ == ny )
        return fineGridOperator_.getL(position,sx,sy,nx,ny);
    const Quadruple index( sx, sy, nx, ny );
        std::cout<<"living"<<std::endl;
    std::map<Quadruple,NumericArray>::const_iterator result = 
                                                    data_[position].find(index);
        std::cout<<"living"<<std::endl;
    if ( result != data_[position].end() )
        return result->second;
                                                
    PositionArray jX;
    PositionArray jY;
    NumericArray operatorL;
        std::cout<<"living"<<std::endl;
    computeGalerkin(
        operatorL,jX,jY,
        position,sx,sy,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
        std::cout<<"living"<<std::endl;
    data_[position].insert(std::make_pair(
        index,
        operatorL ) );
    std::cout<<"living"<<std::endl;
    return data_[position][index];
}

}

#endif /*GALERKIN_H_*/
