#ifndef MULTI_DIM_HPP
#define MULTI_DIM_HPP

#include <array>
#include <vector>
#include <cassert>
#include <stdexcept>

namespace multi_dim{
    using std::array;
    using std::vector;
    using std::range_error;

    template<typename ElementType_, std::size_t SIZE1_, std::size_t SIZE2_>
    class Array2D
    {
      
        private:
            typedef array<ElementType_, SIZE1_ * SIZE2_> Container_;

        public:
            typedef typename Container_::value_type      value_type;
            typedef typename Container_::reference       reference;
            typedef typename Container_::const_reference const_reference;
            typedef typename Container_::size_type       size_type;

        public:
            reference at(const size_type i, const size_type j)
            {
                rangeCheck_(i, j);
                return buffer_[i * SIZE2_ + j];
            }

            const_reference at(const size_type i, const size_type j) const
            {
                rangeCheck_(i, j);
                return buffer_[i * SIZE2_ + j];
            }

            inline reference operator()(const size_type i, const size_type j)
            {
                return buffer_[i * SIZE2_ + j];
            }

            inline const_reference operator()(const size_type i, const size_type j) const
            {
                return buffer_[i * SIZE2_ + j];
            }

            size_type size1() const
            {
                return SIZE1_;
            }
            size_type size2() const
            {
                return SIZE2_;
            }

          private:
            void rangeCheck_(const size_type i, const size_type j) const
            {
                if((i >= SIZE1_) || (j >= SIZE2_))
                    throw range_error();
            }
            
          private:
            Container_ buffer_;
    };

    template<typename ElementType_>
    class Matrix
    {
        private:
            typedef std::vector<ElementType_>   Container_;
        public:
            typedef typename Container_::value_type              value_type;
            typedef typename Container_::reference               reference;
            typedef typename Container_::const_reference         const_reference;
            typedef typename Container_::size_type               size_type;
        public:
            Matrix(const size_type nrows, const size_type ncols)
              : buffer_(nrows * ncols), size1_(nrows), size2_(ncols)
            {}

        #if MATRIX_RESIZABLE
      
            Matrix()
                : buffer_(), size1_(0), size2_(0)
                {}

            void clear()
            {
                size1_ = 0;
                size2_ = 0;
                buffer_.clear();
            }

            void resize(const size_type nrows, const size_type ncols)
            {
                size1_ = nrows;
                size2_ = ncols;
                buffer_.resize(nrows * ncols);
            }
        #endif // MATRIX_RESIZABLE

        reference at(const size_type i, const size_type j)
        {
            rangeCheck_(i, j);
            return buffer_[i * size2_ + j];
        }
        
        const_reference at(const size_type i, const size_type j) const
        {
            rangeCheck_(i, j);
            return buffer_[i * size2_ + j];
        }

        inline reference operator()(const size_type i, const size_type j)
        {
            return buffer_[i * size2_ + j];
        }
        
        inline const_reference operator()(const size_type i, const size_type j) const
        {
            return buffer_[i * size2_ + j];
        }

        size_type size1() const
        {
            return size1_;
        }
      
      size_type size2() const
      {
          return size2_;
      }

      private:
        void rangeCheck_(const size_type i, const size_type j) const
        {
            if((i >= size1_) || (j >= size2_))
                throw range_error();
        }
        
      private:
          Container_ buffer_;
          size_type size1_, size2_;
    };

  }

#endif
