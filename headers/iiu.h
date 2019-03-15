/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */

#ifndef IIU_H_
#define IIU_H_

#include "common.h"
#include "codecs.h"
#include "bitpacking.h"
#include "util.h"

namespace FastPForLib {

/**
 * This implements as best as possible the PFor scheme
 * from  Zukowski et al., Super-Scalar RAM-CPU Cache Compression.
 *
 *  Implemented by D. Lemire
 *
 * In a multithreaded context, you may need one PFor per thread.
 *
 * Small differences:
 *  1. We don't write the exception section is reverse order.
 *     It is written in forward order.
 *
 *  2.  Obviously, the code is specific to 32-bit integers whereas
 *      the original description allowed for different data types.
 *
 *  3. Because we assume pfor delta, we don't compute the base (that is
 *     the frame). Correspondingly, we don't sort a sample. Instead,
 *     we use a fast approach to identify the best number of bits
 *     based on the computation of the integer logarithm.
 *
 *  4. Though it is not clear how the sample is taken in the original
 *     paper, we consider a consecutive sample of up to 64K samples.
 *
 */
class IIU : public IntegerCODEC {
public:
  enum {
    BlockSizeInUnitsOfPackSize = 4,
    PACKSIZE = 32,
    BlockSize = BlockSizeInUnitsOfPackSize * PACKSIZE,
    blocksizeinbits = 7 // constexprbits(BlockSize)
  };
  // these are reusable buffers
  std::vector<uint32_t> codedcopy;
  std::vector<uint32_t> miss;
  typedef uint32_t
      DATATYPE; // this is so that our code looks more like the original paper

  IIU() : codedcopy(BlockSize), miss(BlockSize) {}
  // We use static b bit
  static uint32_t determineBestBase(const DATATYPE *in, size_t size) {
    if (size == 0)
      return 0;
    //const size_t defaultsamplesize = 64 * 1024;
    //// the original paper describes sorting
    //// a sample, but this only makes sense if you
    //// are coding a frame of reference.
    //size_t samplesize = size > defaultsamplesize ? defaultsamplesize : size;
    //uint32_t freqs[33];
    //for (uint32_t k = 0; k <= 32; ++k)
    //  freqs[k] = 0;
    //// we choose the sample to be consecutive
    //uint32_t rstart =
    //    size > samplesize
    //        ? (rand() % (static_cast<uint32_t>(size - samplesize)))
    //        : 0U;
    //for (uint32_t k = rstart; k < rstart + samplesize; ++k) {
    //  freqs[asmbits(in[k])]++;
    //}
    uint32_t bestb = 7;	// 3, 7, 15
    //uint32_t numberofexceptions = 0;
    //double Erate = 0;
    //double bestcost = 32;
    //for (uint32_t b = bestb - 1; b < 32; --b) {
    //  numberofexceptions += freqs[b + 1];
    //  Erate = static_cast<double>(numberofexceptions) /
    //          static_cast<double>(samplesize);
    //  /**
    //  * though this is not explicit in the original paper, you
    //  * need to somehow compensate for compulsory exceptions
    //  * when the chosen number of bits is small.
    //  *
    //  * We use their formula (3.1.5) to estimate actual number
    //  * of total exceptions, including compulsory exceptions.
    //  */
    //  if (numberofexceptions > 0) {
    //    double altErate = (Erate * 128 - 1) / (Erate * (1U << b));
    //    if (altErate > Erate)
    //      Erate = altErate;
    //  }
    //  const double thiscost = b + Erate * 32;
    //  if (thiscost <= bestcost) {
    //    bestcost = thiscost;
    //    bestb = b;
    //  }
    //}
    return bestb;
  }

  // returns location of first exception or BlockSize if there is none
  uint32_t compressblockPFOR(const DATATYPE *__restrict__ in,
                             uint32_t *__restrict__ outputbegin,
                             const uint32_t b,
                             //DATATYPE *__restrict__ &exceptions) {
                             DATATYPE *__restrict__ &exceptions, uint32_t *__restrict__ skiplist) {
    if (b == 32) {
      for (size_t k = 0; k < BlockSize; ++k)
        *(outputbegin++) = *(in++);
      return BlockSize;
    }
    size_t exceptcounter = 0;
    const uint32_t maxgap = 1U << b;
		uint32_t exceptcounter1[2] = {0, 0}; // for skiplist
		uint32_t z = 0;
    {
      std::vector<uint32_t>::iterator cci = codedcopy.begin();
			//std::cout << "[arcj] in: " << std::endl;
      for (uint32_t k = 0; k < BlockSize; ++k, ++cci) {
				//std::cout << "(" << k << "/" << in[k] << "/" << (in[k] >= maxgap) << ") ";
        miss[exceptcounter] = k; // the position of an exception value
        exceptcounter += (in[k] >= maxgap);
				if (k % 64 == 63) {	exceptcounter1[z++] = exceptcounter; }
      }
			*(skiplist + 1) = exceptcounter1[0];
			*(skiplist + 5) = exceptcounter1[1] - exceptcounter1[0];
    }
    if (exceptcounter == 0) {
			//packblock(in, outputbegin, b);
      packblock(in, outputbegin, b, &miss[0], exceptcounter);
      return BlockSize;
    }
    codedcopy.assign(in, in + BlockSize);
    uint32_t firstexcept = miss[0];
    //uint32_t prev = 0;
		size_t tmp = 0;
    *(exceptions++) = codedcopy[firstexcept]; // store the first exception value
		//codedcopy[firstexcept] = firstexcept + 128;
		codedcopy[firstexcept] = tmp + 128; // set 8-th bit (exception bit)
		tmp++;
    //prev = firstexcept;
    //if (maxgap < BlockSize) {
    //  for (uint32_t i = 1; i < exceptcounter; ++i) {
    //    uint32_t cur = miss[i];
    //    // they don't include this part, but it is required:
    //    while (cur > maxgap + prev) {
    //      // compulsory exception
    //      uint32_t compulcur = prev + maxgap;
    //      *(exceptions++) = codedcopy[compulcur];
    //      codedcopy[prev] = maxgap - 1;
    //      prev = compulcur;
    //    }
    //    *(exceptions++) = codedcopy[cur];
    //    codedcopy[prev] = cur - prev - 1;
    //    prev = cur;
    //  }
    //} else {
      for (uint32_t i = 1; i < exceptcounter; ++i) {
        uint32_t cur = miss[i];
        *(exceptions++) = codedcopy[cur];
				codedcopy[cur] = tmp + 128;
				tmp++;
        //codedcopy[prev] = cur - prev - 1; // arcj: the offset of next exception from this value
				//prev = cur;
      }
    //}
		//	std::cout << "[arcj] codedcopy: " << std::endl;
		//for (uint32_t z = 0; z < BlockSize; z++) {
		//	std::cout << codedcopy[z] << "/";
		//	if(z % 64 == 63) { std::cout << std::endl; }
		//}
		//std::cout << std::endl;
    //packblock(&codedcopy[0], outputbegin, b); // arcj: org
    packblock(&codedcopy[0], outputbegin, b, &miss[0], exceptcounter);
    return firstexcept; // arcj: the position (offset) of first exception value
  }

  // void packblock(const uint32_t *source, uint32_t *out, const uint32_t bit) {
  void packblock(const uint32_t *source, uint32_t *out, const uint32_t bit /* b bit */, const uint32_t *exception, const size_t exceptcounter) {
    for (uint32_t j = 0; j != BlockSize; j += PACKSIZE) { // per 32 elements
      fastpack_iiu(source + j, out, bit, exception, exceptcounter, j);
      //out += bit;
      out += (bit + 1);
    }
  }

  void unpackblock(const uint32_t *source, uint32_t *out, const uint32_t bit) {
    for (uint32_t j = 0; j != BlockSize; j += PACKSIZE) {
      fastunpack_iiu(source, out + j, bit);
      //source += bit;
      source += (bit + 1);
    }
  }

  void encodeArray(const uint32_t *in, const size_t len, uint32_t *out,
                   //size_t &nvalue) {
                   size_t &nvalue, uint32_t *skiplist) {
    *out++ = static_cast<uint32_t>(len);
#ifndef NDEBUG
    const uint32_t *const finalin(in + len);
#endif
    const uint32_t maxsize = (1U << (32 - blocksizeinbits - 1)); // blocksizeinbits = 7, 128 entries per block
    size_t totalnvalue(1);
    // for (size_t i = 0; i < len; i += maxsize)
    for (size_t j = 0; j < (len + maxsize - 1U) / maxsize; ++j) {
      size_t i = j << (32 - blocksizeinbits - 1); // 32 - 7 - 1 = 24
      size_t l = maxsize;
      if (i + maxsize > len) {
        l = len - i;
        assert(l <= maxsize);
      }
      size_t thisnvalue = nvalue - totalnvalue;
      assert(in + i + l <= finalin);
      //__encodeArray(&in[i], l, out, thisnvalue);
      __encodeArray(&in[i], l, out, thisnvalue, skiplist);
      totalnvalue += thisnvalue;
      assert(totalnvalue <= nvalue);
      out += thisnvalue;
    }
    nvalue = totalnvalue;
  }
  const uint32_t *decodeArray(const uint32_t *in, const size_t len,
                              //uint32_t *out, size_t &nvalue) {
                              uint32_t *out, size_t &nvalue, uint32_t *skiplist) {
    nvalue = *in++;
    if (nvalue == 0)
      return in;
#ifndef NDEBUG
    const uint32_t *const initin = in;
#endif
    const uint32_t *const finalin = in + len;
    size_t totalnvalue(0);
    while (totalnvalue < nvalue) {
      size_t thisnvalue = nvalue - totalnvalue;
#ifndef NDEBUG
      const uint32_t *const befin(in);
#endif
      assert(finalin <= len + in);
      //in = __decodeArray(in, finalin - in, out, thisnvalue);
      in = __decodeArray(in, finalin - in, out, thisnvalue, skiplist);
      assert(in > befin);
      assert(in <= finalin);
      out += thisnvalue;
      totalnvalue += thisnvalue;
      assert(totalnvalue <= nvalue);
    }
    assert(in <= len + initin);
    assert(in <= finalin);
    nvalue = totalnvalue;
    return in;
  }

  void __encodeArray(const uint32_t *in, const size_t len, uint32_t *out,
                     //size_t &nvalue) {
                     size_t &nvalue, uint32_t *skiplist) {
    checkifdivisibleby(len, BlockSize);
    const uint32_t *const initout(out);
    std::vector<DATATYPE> exceptions; // uint32_t
    exceptions.resize(len);
    DATATYPE *__restrict__ i = &exceptions[0]; // 1 exception vector
    const uint32_t b = determineBestBase(in, len); 
    *out++ = static_cast<uint32_t>(len);
    *out++ = b;
		
		// arcj: mod
    for (size_t k = 0; k < len / BlockSize; ++k) { // BlockSize = 128
			*(skiplist + 2) = size_t(out) >> 32;
			*(skiplist + 3) = size_t(out) & 0xFFFFFFFF;
		  *(skiplist + 6) = size_t(out + (BlockSize * (b + 1)) / 32 / 2) >> 32;
		  *(skiplist + 7) = size_t(out + (BlockSize * (b + 1)) / 32 / 2) & 0xFFFFFFFF;
      //uint32_t firstexcept = compressblockPFOR(in, out, b, i, skiplist);
			DATATYPE *tmp = i;
			compressblockPFOR(in, out, b, i, skiplist);
			// check values in skiplist
			//std::cout << "[arcj] skiplist1: " << *(skiplist + 0) << "/" << *(skiplist + 1) << "/" << *(skiplist + 2) << "/" << *(skiplist + 3) << "/" << std::endl;
			//std::cout << "[arcj] skiplist2: " << *(skiplist + 4) << "/" << *(skiplist + 5) << "/" << *(skiplist + 6) << "/" << *(skiplist + 7) << "/" << std::endl;
      out += (BlockSize * (b + 1)) / 32; 
      in += BlockSize;
			for (uint32_t t = 0; t < (*(skiplist + 1) + *(skiplist + 5)); ++t) {
		  	*out++ = *(tmp + t);
			}
			skiplist += 8;
    }
    nvalue = out - initout;

		//// arcj: org
    //for (size_t k = 0; k < len / BlockSize; ++k) { // arcj: BlockSize = 128
    //  uint32_t *const headerout(out);
    //  ++out;
    //  //uint32_t firstexcept = compressblockPFOR(in, out, b, i); // arcj: org
		//	//*(skiplist + 2) = size_t(out) >> 32;
		//	//*(skiplist + 3) = size_t(out) & 0xFFFFFFFF;
		//	//*(skiplist + 6) = size_t(out + (BlockSize * (b + 1)) / 32 / 2) >> 32;
		//	//*(skiplist + 7) = size_t(out + (BlockSize * (b + 1)) / 32 / 2) & 0xFFFFFFFF;
		//	//std::cout << "[arcj] out1 " << out << "/" << out + (BlockSize * (b + 1)) / 32 / 2 << std::endl;
		//	//std::cout << "[arcj] encode1 " << i << std::endl;
    //  uint32_t firstexcept = compressblockPFOR(in, out, b, i, skiplist);
    //  out += (BlockSize * (b + 1)) / 32; 
		//	// out += (BlockSize * b) / 32; // arcj: org
		//	std::cout << "[arcj] encode1 " << *(skiplist + 5) + *(skiplist + 1) << std::endl;
    //  in += BlockSize;
		//	skiplist += 8;
		//	// arcj: we don't need this part
    //  const uint32_t bitsforfirstexcept = blocksizeinbits;
    //  const uint32_t firstexceptmask = (1U << blocksizeinbits) - 1;
    //  const uint32_t exceptindex = static_cast<uint32_t>(i - &exceptions[0]);
		//	std::cout << "[arcj] encode2 " << i  << "/" << exceptindex << std::endl;
    //  *headerout =
    //      (firstexcept & firstexceptmask) | (exceptindex << bitsforfirstexcept);
    //}
    //const size_t howmanyexcept = i - &exceptions[0];
		////std::cout << "[arcj] encode3 " <<howmanyexcept << "/" << i << "/" << &exceptions[0] << std::endl;
    //for (uint32_t t = 0; t < howmanyexcept; ++t)
    //  *out++ = exceptions[t];
    //nvalue = out - initout;
  }

#ifndef NDEBUG
  const uint32_t *__decodeArray(const uint32_t *in, const size_t len,
#else
  const uint32_t *__decodeArray(const uint32_t *in, const size_t,
#endif
                                //uint32_t *out, size_t &nvalue) {
                                uint32_t *out, size_t &nvalue, uint32_t *skiplist) {
#ifndef NDEBUG
    const uint32_t *const initin(in);
#endif
    nvalue = *in++;
    checkifdivisibleby(nvalue, BlockSize);
    const uint32_t b = *in++;
    const uint32_t *except = in;
    //const DATATYPE *__restrict__ except =
        //in + nvalue * b / 32 + nvalue / BlockSize; // arcj: org 
    //    in + nvalue * (b + 1) / 32 + nvalue / BlockSize;
    //const uint32_t bitsforfirstexcept = blocksizeinbits;
    //const uint32_t firstexceptmask = (1U << blocksizeinbits) - 1;
    //const DATATYPE *endexceptpointer;
    const uint32_t *endexceptpointer = in;
    //const DATATYPE *const initexcept(except);
    uint32_t arc_prev = 0;

		// arcj: mod
    for (size_t k = 0; k < nvalue / BlockSize; ++k) {
      const uint32_t firstexcept = 0; // we don't need it
      const uint32_t exceptindex = 0;
			// calculate exception address
			except += ((BlockSize * (b + 1)) / 32);
			endexceptpointer += ((BlockSize * (b + 1)) / 32) + *(skiplist + 1) + *(skiplist + 5);
      uncompressblockPFOR(in, out, b, except, /*we don't need following args*/endexceptpointer, firstexcept, exceptindex, arc_prev, skiplist);
      in += (BlockSize * (b + 1)) / 32 + *(skiplist + 1) + *(skiplist + 5);
      out += BlockSize;
			skiplist += 8;
		}

		// arcj: org
    //for (size_t k = 0; k < nvalue / BlockSize; ++k) {
    //  const uint32_t *const headerin(in);
    //  ++in;
    //  const uint32_t firstexcept = *headerin & firstexceptmask;
    //  const uint32_t exceptindex = *headerin >> bitsforfirstexcept;
    //  endexceptpointer = initexcept + exceptindex;
		//	//std::cout << "[arcj] decode " << firstexcept << "/" << exceptindex << "/" << endexceptpointer << std::endl;
    //  //uncompressblockPFOR(in, out, b, except, endexceptpointer, firstexcept); // arcj: org
    //  uncompressblockPFOR(in, out, b, except, endexceptpointer, firstexcept, exceptindex, arc_prev, skiplist);
		//	arc_prev = exceptindex;
    //  //in += (BlockSize * b) / 32; // arcj: org
    //  in += (BlockSize * (b + 1)) / 32;
    //  out += BlockSize;
		//	skiplist += 8;
    //}
    assert(initin + len >= in);
    assert(initin + len >= endexceptpointer);
    return endexceptpointer;
  }
  void uncompressblockPFOR(
      const uint32_t
          *__restrict__ inputbegin, // points to the first packed word
      DATATYPE *__restrict__ outputbegin,
      const uint32_t b,
      //const DATATYPE *__restrict__
      const uint32_t *
          &i, // i points to value of the first exception
      const uint32_t *end_exception,
      //const DATATYPE *__restrict__ end_exception,
      size_t next_exception, // points to the position of the first exception
      const uint32_t exceptindex,
			const uint32_t prev,
			uint32_t *__restrict__ skiplist
      ) {
    unpackblock(inputbegin, reinterpret_cast<uint32_t *>(outputbegin),
                b); /* bit-unpack the values */
		for (uint32_t z = 0; z < 128; ++z) {
			if (outputbegin[z] >= 128) {
				outputbegin[z] = *(i++);
			} 
		}
		// arcj: org
    //for (size_t cur = next_exception; i != end_exception;
    //     cur = next_exception) {
    //  next_exception = cur + static_cast<size_t>(outputbegin[cur]) + 1;
    //  outputbegin[cur] = *(i++);
    //}
  }

  virtual std::string name() const {
    std::ostringstream convert;
    convert << "IIU";
    return convert.str();
  }
};

} // namespace FastPFor

#endif /* IIU_H_ */
