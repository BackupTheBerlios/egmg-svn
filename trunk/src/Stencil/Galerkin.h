#ifndef GALERKIN_H_
#define GALERKIN_H_

#include<map>
#include<vector>
#include<valarray>
#include "Stencil.h"
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

class Galerkin : public mg::Stencil
{
private:
	class Quadruple
	{
	private:
		const size_t i_,j_,k_,l_;
	public:
		Quadruple(const size_t i,const size_t j,const size_t k,const size_t l)
			:i_(i),j_(j),k_(k),l_(l) {}
		/*Quadruple(const Quadruple& rhs)
			:i_(rhs.i_),j_(rhs.j_),k_(rhs.k_),l_(rhs.l_) {}
		~Quadruple() {}
		Quadruple& operator=(const Quadruple& rhs)
		{
			if (this == &rhs)
				return *this;
			i_=rhs.i_;
			j_=rhs.j_;
			k_=rhs.k_;
			l_=rhs.l_;
			return *this;
		}*/
		bool operator==(const Quadruple& rhs) const
		{
			if (this == &rhs)
				return true;
			return i_==rhs.i_ && j_==rhs.j_	&& k_== rhs.k_ && l_==rhs.l_;
		}
		bool operator!=(const Quadruple& rhs) const
		{
			if (this == &rhs)
				return false;
			return !(*this==rhs);
		}
		bool operator>(const Quadruple& rhs) const
		{
			return !(*this<rhs);
		}
		bool operator<(const Quadruple& rhs) const
		{
			
			return (i_<rhs.i_||i_==rhs.i_&&j_<rhs.j_
							 ||i_==rhs.i_&&j_==rhs.j_&&k_<rhs.K_
							 ||i_==rhs.i_&&j_==rhs.j_
							 			 &&k_==rhs.k_&&l_<rhs.l_);
		}
	};
	std::vector<std::map<Quadruple,std::valarray<precision> > > data;
	std::vector<std::valarray<int> > J_x;
	std::vector<std::valarray<int> > J_y;
	const Stencil& stencil;
	size_t size_;
	std::vector<const Prolongation*> prolongs;
	std::vector<const Restriction*> restricts;
	std::vector<std::valarray<int> > init_J_x(const Stencil& sten);
	std::vector<std::valarray<int> > init_J_y(const Stencil& sten);
	void update_size();
	void update_J_x_J_y();
public:
	Galerkin(const Stencil& fine_grid_operator)
		: J_x(init_J_x(fine_grid_operator)),J_y(init_J_y(fine_grid_operator))
			,stencil(fine_grid_operator),size_(fine_grid_operator.size()) {}
	virtual ~Galerkin() {}
	virtual precision apply_c(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_w(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_nw(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_n(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_ne(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_e(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_se(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_s(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision apply_sw(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const;
	virtual precision get_center_c(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_w(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_nw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_n(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_ne(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_e(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_se(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_s(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual precision get_center_sw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_c(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		if (data[c].find(Quadruple(i,j,Nx,Ny)) == data[c].end())
				data[c].insert(Quadruple(i,j,Nx,Ny),compute_L(c,i,j,Nx,Ny);
		return data[c][Quadruple(i,j,Nx,Ny)];
	}
	virtual const std::valarray<precision>& get_L_w(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_nw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_n(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_ne(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_e(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_se(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_s(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	virtual const std::valarray<precision>& get_L_sw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const;
	inline const std::valarray<int>& get_J_x(const pos p =c) const
	{
		return J_x[p];
	}
	inline const std::valarray<int>& get_J_y(const pos p =c) const
	{
		return J_y[p];
	}
	void push_prolongation(const Prolongation* prolongation)
	{
		prolongs.push_back(prolongation);
		update_size();
		update_J_x_J_y();
	}
	void pop_prolongation()
	{
		prolongs.pop_back();
		update_size();
		update_J_x_J_y();
	}
	void push_restriction(const Restriction* restriction)
	{
		restricts.push_back(restriction);
		update_size();
		update_J_x_J_y();
	}
	void pop_restriction()
	{
		restricts.pop_back();
		update_size();
		update_J_x_J_y();
	}
	inline size_t size() const
	{
		return size_;
	}
	inline bool is_constant() const
	{
		return stencil.is_constant();
	}
};

}

#endif /*GALERKIN_H_*/
