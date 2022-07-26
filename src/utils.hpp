/*
 Crown Copyright 2012 AWE.

 This file is part of CloverLeaf.

 CloverLeaf is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the
 Free Software Foundation, either version 3 of the License, or (at your option)
 any later version.

 CloverLeaf is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details.

 You should have received a copy of the GNU General Public License along with
 CloverLeaf. If not, see http://www.gnu.org/licenses/.
 */

#ifndef UTILS_HPP
#define UTILS_HPP

#include "dpl_shim.h"

#include <iostream>
#include <utility>
#include <cassert>
#include <vector>

namespace clover {


	struct Range1d {
		const size_t from, to;
		template<typename A, typename B>
		Range1d(A from, B to) : from(from), to(to) {
			assert(from < to);
			assert(to - from > 0);
		}
		friend std::ostream &operator<<(std::ostream &os, const Range1d &d) {
			os << "Range1d{"
			   << " X[" << d.from << "->" << d.to << " (" << (d.to - d.from) << ")]"
			   << "}";
			return os;
		}
	};

	struct Range2d {
		const size_t fromX, toX;
		const size_t fromY, toY;
		const size_t sizeX, sizeY;
		template<typename A, typename B, typename C, typename D>
		Range2d(A fromX, B fromY, C toX, D toY) :
				fromX(fromX), toX(toX), fromY(fromY), toY(toY),
				sizeX(toX - fromX), sizeY(toY - fromY) {
			assert(fromX < toX);
			assert(fromY < toY);
			assert(sizeX != 0);
			assert(sizeY != 0);
		}
		friend std::ostream &operator<<(std::ostream &os, const Range2d &d) {
			os << "Range2d{"
			   << " X[" << d.fromX << "->" << d.toX << " (" << d.sizeX << ")]"
			   << " Y[" << d.fromY << "->" << d.toY << " (" << d.sizeY << ")]"
			   << "}";
			return os;
		}
	};

    template <typename N>
    class range {
    public:
        class iterator {
            friend class range;
        public:

            using difference_type = typename std::make_signed_t<N>;
            using value_type = N;
            using pointer = const N*;
            using reference = N;
            using iterator_category = std::random_access_iterator_tag;

            // XXX This is not part of the iterator spec, it gets picked up by oneDPL if enabled.
            // Without this, the DPL SYCL backend collects the iterator data on the host and copies to the device.
            // This type is unused for any other STL impl.
            using is_passed_directly = std::true_type;

            reference operator *() const { return i_; }
            iterator &operator ++() { ++i_; return *this; }
            iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }

            iterator &operator --() { --i_; return *this; }
            iterator operator --(int) { iterator copy(*this); --i_; return copy; }

            iterator &operator +=(N by) { i_+=by; return *this; }

            value_type operator[](const difference_type &i) const { return i_ + i; }

            difference_type operator-(const iterator &it) const { return i_ - it.i_; }
            iterator operator+(const value_type v) const { return iterator(i_ + v); }

            bool operator ==(const iterator &other) const { return i_ == other.i_; }
            bool operator !=(const iterator &other) const { return i_ != other.i_; }
            bool operator < (const iterator &other) const { return i_ < other.i_; }

        protected:
            explicit iterator(N start) : i_ (start) {}

        private:
            N i_;
        };

        iterator begin() const { return begin_; }
        iterator end() const { return end_; }
        range(N begin, N end) : begin_(begin), end_(end) {}
    private:
        iterator begin_;
        iterator end_;
    };

    // seq | par | par_unseq | unseq
//    static auto EXEC_POLICY = exe_policy;

    template<typename F>
	static void par_ranged1(const Range1d &r, const F &functor) {
		auto groups = range<int>(r.from, r.to);
		std::for_each(EXEC_POLICY,  groups.begin(), groups.end(), [functor](int i) {
			functor(i);
		});
		//for (size_t i = r.from; i < r.to; i++) {
		//	functor(i);
		//}
	}

	template<typename F>
	static void par_ranged2(const Range2d &r, const F &functor) {
        auto xy = range<int>(0, r.sizeX * r.sizeY);
        std::for_each(EXEC_POLICY, xy.begin(), xy.end(), [=](int v) {
            const auto x = r.fromX + (v % r.sizeX);
            const auto y = r.fromY + (v / r.sizeX);
            functor(x, y);
        });
        //for (size_t j = r.fromY; j < r.toY; j++) {
        //    for (size_t i = r.fromX; i < r.toX; i++) {
        //        functor(i, j);
        //    }
        //}
	}

#ifdef USE_VECTOR

    template<typename T>
    struct Buffer1D {

        const size_t size;
        std::vector<T> data;

        explicit Buffer1D(size_t size) : size(size), data(size) {}
        T operator[](size_t i) const { return data[i]; }
        T &operator[](size_t i) { return data[i]; }
        T *actual() { return data.data(); }

        friend std::ostream &operator<<(std::ostream &os, const Buffer1D<T> &buffer) {
            os << "Buffer1D(size: " << buffer.size << ")";
            return os;
        }

    };

    template<typename T>
    struct Buffer2D {

        const size_t sizeX, sizeY;
        std::vector<T> data;

        Buffer2D(size_t sizeX, size_t sizeY) : sizeX(sizeX), sizeY(sizeY), data(sizeX * sizeY) {}
        T &operator()(size_t i, size_t j) { return data[i + j * sizeX]; }
        T const &operator()(size_t i, size_t j) const { return data[i + j * sizeX]; }
        T *actual() { return data.data(); }

        friend std::ostream &operator<<(std::ostream &os, const Buffer2D<T> &buffer) {
            os << "Buffer2D(sizeX: " << buffer.sizeX << " sizeY: " << buffer.sizeY << ")";
            return os;
        }

    };

#else


    template<typename T>
    struct Buffer1D {

        size_t size;
        T *data;

        explicit Buffer1D(size_t size) : size(size), data(alloc_raw<T>( size)) {}
//        Buffer1D(const Buffer1D &that) : Buffer1D(that.size) {
//            // XXX parallel copy should work on GPUs, but we get illegal memory access in NVHPC stdpar. Just copy on host for now, this is only needed for decomposition (it's not timed).
//#ifdef SERIAL_COPY_CTOR
//            std::copy(that.data, that.data + size, data);
//#else
//            par_ranged1({0, size}, [=, this](auto i){ data[i] = that.data[i]; });
//#endif
//        }
//        Buffer1D &operator=(const Buffer1D &other) = delete;
//        T operator[](size_t i) const { return data[i]; }
        T &operator[](size_t i)  const { return data[i]; }
        T *actual() { return data; }
        //virtual ~Buffer1D() { std::free(data); }

        friend std::ostream &operator<<(std::ostream &os, const Buffer1D<T> &buffer) {
            os << "Buffer1D(size: " << buffer.size << ")";
            return os;
        }

    };

    template<typename T>
    struct Buffer2D {

        size_t sizeX, sizeY;
        T *data;

        Buffer2D(size_t sizeX, size_t sizeY) : sizeX(sizeX), sizeY(sizeY),  data(alloc_raw<T>(sizeX * sizeY)) {}
//        Buffer2D(const Buffer2D &that) : Buffer2D(that.sizeX, that.sizeY) {
//            // XXX parallel copy should work on GPUs, but we get illegal memory access in NVHPC stdpar. Just copy on host for now, this is only needed for decomposition (it's not timed).
//#ifdef SERIAL_COPY_CTOR
//            std::copy(that.data, that.data + that.sizeX * that.sizeY, data);
//#else
//            par_ranged2({0, 0, sizeX, sizeY}, [=, this](auto i, auto j){ data[i + j * sizeX] = that.data[i + j * sizeX]; });
//#endif
//        }
//        Buffer2D &operator=(const Buffer2D &other) = delete;
        T &operator()(size_t i, size_t j) const { return data[i + j * sizeX]; }
//        T &operator()(size_t i, size_t j) const { return data[i + j * sizeX]; }
        T *actual() { return data; }
        //virtual ~Buffer2D() { std::free(data); }

        friend std::ostream &operator<<(std::ostream &os, const Buffer2D<T> &buffer) {
            os << "Buffer2D(sizeX: " << buffer.sizeX << " sizeY: " << buffer.sizeY << ")";
            return os;
        }
    };

#endif


//    template<typename T>
//    struct Buffer1DX {
//
//        size_t size;
//        T *data;
//
//        explicit Buffer1DX(size_t size) : size(size), data(alloc_raw<T>( size)) {}
//
//
//
//        Buffer1DX &operator=(const Buffer1DX &other) = delete;
//        T operator[](size_t i) const { return data[i]; }
//        T &operator[](size_t i) { return data[i]; }
//        T *actual() { return data; }
//        //virtual ~Buffer1D() { std::free(data); }
//
//        friend std::ostream &operator<<(std::ostream &os, const Buffer1D<T> &buffer) {
//            os << "Buffer1D(size: " << buffer.size << ")";
//            return os;
//        }
//
//    };
//
//    static_assert(sycl::is_device_copyable_v<Buffer1DX<float>>);


}


using namespace clover;

#endif //UTILS_HPP
