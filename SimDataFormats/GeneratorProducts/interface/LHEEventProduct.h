#ifndef SimDataFormats_GeneratorProducts_LHEEventProduct_h
#define SimDataFormats_GeneratorProducts_LHEEventProduct_h

#include <memory>
#include <vector>
#include <string>

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/WeightsInfo.h"

class LHEEventProduct {
    public:
	typedef gen::PdfInfo PDF;
	typedef gen::WeightsInfo WGT;

	typedef std::vector<std::string>::const_iterator
						comments_const_iterator;
	typedef std::vector<std::string>::size_type size_type;

	LHEEventProduct() {}
	LHEEventProduct(const lhef::HEPEUP &hepeup) : hepeup_(hepeup) {}
	~LHEEventProduct() {}

	void setPDF(const PDF &pdf) { pdf_.reset(new PDF(pdf)); }
	void addWeight(const WGT& wgt) { 	  
	  weights_.push_back(wgt);
	}
	void addComment(const std::string &line) { comments_.push_back(line); }

	const std::vector<WGT>& weights() const { return weights_; }

	const lhef::HEPEUP &hepeup() const { return hepeup_; }
	const PDF *pdf() const { return pdf_.get(); }

	size_type comments_size() const { return comments_.size(); }
	comments_const_iterator comments_begin() const { return comments_.begin(); }
	comments_const_iterator comments_end() const { return comments_.end(); }

	class const_iterator {
	    public:
		typedef std::forward_iterator_tag	iterator_category;
		typedef std::string			value_type;
		typedef std::ptrdiff_t			difference_type;
		typedef std::string			*pointer;
		typedef std::string			&reference;

		const_iterator() : line(npos) {}
		~const_iterator() {}

		inline bool operator == (const const_iterator &other) const
		{ return line == other.line; }
		inline bool operator != (const const_iterator &other) const
		{ return !operator == (other); }

		inline const_iterator &operator ++ ()
		{ next(); return *this; }
		inline const_iterator operator ++ (int dummy)
		{ const_iterator orig = *this; next(); return orig; }

		const std::string &operator * () const { return tmp; }
		const std::string *operator -> () const { return &tmp; }

	    private:
		friend class LHEEventProduct;

		void next();

		const LHEEventProduct	*event;
		unsigned int		line;
		std::string		tmp;

		static const unsigned int npos = 99999;
	};

	const_iterator begin() const;
	inline const_iterator end() const { return const_iterator(); }

    private:
	lhef::HEPEUP			hepeup_;
	std::vector<std::string>	comments_;
	std::auto_ptr<PDF>		pdf_;
	std::vector<WGT>                weights_;
};

#endif // GeneratorEvent_LHEInterface_LHEEventProduct_h
