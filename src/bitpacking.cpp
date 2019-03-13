#include "bitpacking.h"
#include <cstdint>
#include <type_traits>

namespace {

template<uint8_t DELTA, uint8_t SHR> typename std::enable_if<(DELTA + SHR) < 32>::type
  unpack_single_out(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {

  *out = ((*in) >> SHR) % (1 << DELTA);
}

template<uint8_t DELTA, uint8_t SHR>
  typename std::enable_if<(DELTA + SHR) >= 32>::type
  unpack_single_out(const uint32_t *__restrict__ & in, uint32_t *__restrict__ out) {

  *out = (*in) >> SHR;
  ++in;

  static const uint8_t NEXT_SHR = SHR + DELTA - 32;
  *out |= ((*in) % (1U << NEXT_SHR)) << (32 - SHR);
}

template<uint16_t DELTA, uint16_t SHL, uint32_t MASK>
  typename std::enable_if<DELTA + SHL >= 32>::type
  pack_single_in(const uint32_t in, uint32_t *__restrict__& out) {

  *out |= in << SHL;
  ++out;

  if (DELTA + SHL > 32) {
    *out = (in  & MASK) >> (32 - SHL);
  }
}

template<uint16_t DELTA, uint16_t SHL, uint32_t MASK>
typename std::enable_if<DELTA + SHL < 32>::type
  pack_single_in(const uint32_t in, uint32_t *__restrict__ out) {

  if (SHL == 0) {
    *out = in  & MASK;
  } else {
    *out |= (in & MASK) << SHL;
  }
}

// arcj: for iiu FIXME
template<uint8_t DELTA, uint8_t SHR> typename std::enable_if<(DELTA + SHR) < 32>::type
  unpack_single_out_iiu(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {

  *out = ((*in) >> SHR) % (1 << DELTA);
	//std::cout << "[arcj] unpack single out1 " << in << " / " << DELTA << " / " << SHR << " / " << DELTA + SHR << std::endl;
}

template<uint8_t DELTA, uint8_t SHR>
  typename std::enable_if<(DELTA + SHR) >= 32>::type
  unpack_single_out_iiu(const uint32_t *__restrict__ & in, uint32_t *__restrict__ out) {
	//std::cout << "[arcj] unpack single out2 " << in << " / " << DELTA << " / " << SHR << " / " << DELTA + SHR << std::endl;

  *out = (*in) >> SHR;
  ++in;

  static const uint8_t NEXT_SHR = SHR + DELTA - 32;
  *out |= ((*in) % (1U << NEXT_SHR)) << (32 - SHR);
}

template<uint16_t DELTA, uint16_t SHL, uint32_t MASK>
  typename std::enable_if<DELTA + SHL >= 32>::type
  pack_single_in_iiu(const uint32_t in, uint32_t *__restrict__& out) {
	//std::cout << "[arcj] pack single in1 " << in << "/ " << DELTA << " / " << SHL << " / " << DELTA + SHL << " / " << MASK << std::endl;
	// case 2
  *out |= in << SHL;
	//std::cout << "[arcj] OUT: " << *out << std::endl;
	//std::cout << "[arcj] out address2-1: " << out << std::endl;
  ++out;
	//std::cout << "[arcj] out address2-2: " << out << std::endl;

	//// arcj: maybe, we don't use it
  //if (DELTA + SHL > 32) {
  //  *out = (in  & MASK) >> (32 - SHL); // arcj: remain bits
  //}
	//std::cout << "[arcj] pack_single_in1 " << in << " / " << *out << std::endl;
}

// arcj: when we can pack in into 32 bit
template<uint16_t DELTA, uint16_t SHL, uint32_t MASK>
typename std::enable_if<DELTA + SHL < 32>::type
  pack_single_in_iiu(const uint32_t in, uint32_t *__restrict__ out) {
	//std::cout << "[arcj] pack single in2 " << in << " / " << DELTA << " / " << SHL << " / " << DELTA + SHL << " / " << MASK << std::endl;

	// case 1
  if (SHL == 0) {
    *out = in  & MASK; // 1111111 (127) FIXME: concate((in & MASK) 1(exception) or 0(normal))
  } else {
    *out |= (in & MASK) << SHL; // FIXME: SHL to SHL + 1bit (8bit)
  }
	//std::cout << "[arcj] pack_single_in2 " << in << " / " << *out << std::endl;
	//std::cout << "[arcj] OUT: " << *out << std::endl;
}


template<uint16_t DELTA, uint16_t OINDEX = 0> struct Unroller {
  static_assert(DELTA < 32, "");

  static void Unpack(const uint32_t *__restrict__ & in, uint32_t *__restrict__ out) {
		//std::cout << "[arcj] Unpack " << DELTA << std::endl;
    unpack_single_out<DELTA, (DELTA * OINDEX) % 32>(in, out + OINDEX);

    Unroller<DELTA, OINDEX + 1>::Unpack(in, out);
  }

  static void Pack(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
    pack_single_in<DELTA, (DELTA * OINDEX) % 32, (1U << DELTA) - 1 >(in[OINDEX], out);

    Unroller<DELTA, OINDEX + 1>::Pack(in, out);
  }

  static void PackNoMask(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
    pack_single_in<DELTA, (DELTA * OINDEX) % 32, uint32_t(-1)>(in[OINDEX], out);

    Unroller<DELTA, OINDEX + 1>::PackNoMask(in, out);
  }

	// arcj: for iiu FIXME
  static void Unpack_iiu(const uint32_t *__restrict__ & in, uint32_t *__restrict__ out) {
		//std::cout << "[arcj] Unpack_iiu " << DELTA << " / " << (DELTA * OINDEX) % 32 << std::endl;
    //unpack_single_out_iiu<DELTA, (DELTA * OINDEX) % 32>(in, out + OINDEX);
    unpack_single_out_iiu<DELTA + 1, ((DELTA + 1) * OINDEX) % 32>(in, out + OINDEX);

    //Unroller<DELTA, OINDEX + 1>::Unpack_iiu(in, out);
    Unroller<DELTA, OINDEX + 1>::Unpack_iiu(in, out);
  }

  static void Pack_iiu(const uint32_t *__restrict__ in, uint32_t *__restrict__ out, const uint32_t *__restrict__ exception, const size_t exceptcounter, const uint32_t count) {
		// arcj: recursive function call
		//std::cout << "[arcj] Pack_iiu1 " << DELTA << " / " << OINDEX << " / " << in[OINDEX] << std::endl;
		//bool exceptiontrue = false;
		//uint32_t pos = 0;
		//for(uint32_t z = 0; z < exceptcounter; ++z) {
		//	if(exception[z] == (OINDEX + count)) {
		//		exceptiontrue = true;
		//		pos = z;
		//	}
		//}
    pack_single_in_iiu<DELTA+1, ((DELTA+1) * OINDEX) % 32, (1U << (DELTA+1)) - 1 >(in[OINDEX], out); // DELTA >> DELTA + 1 

		//if(exceptiontrue) {
		//	//std::cout << "[arcj] Pack_iiu1 exception-" << OINDEX + count << " / "<< exception[pos] << " / " << pos << " >> " << pos + 128 << std::endl;
    //  pack_single_in_iiu<DELTA+1, ((DELTA+1) * OINDEX) % 32, (1U << (DELTA+1)) - 1 >(pos + 128, out); // DELTA >> DELTA + 1 
		//} else {
		//	//std::cout << "[arcj] Pack_iiu1 normal-" << OINDEX + count << " / " << in[OINDEX] << std::endl;
    //  pack_single_in_iiu<DELTA+1, ((DELTA+1) * OINDEX) % 32, (1U << (DELTA+1)) - 1 >(in[OINDEX], out); // DELTA >> DELTA + 1 
		//}
    //pack_single_in_iiu<DELTA, (DELTA * OINDEX) % 32, (1U << DELTA) - 1 >(in[OINDEX], out); // original

    Unroller<DELTA, OINDEX + 1>::Pack_iiu(in, out, exception, exceptcounter, count);
  }
};

// arcj: packing 32th element Done.
template<uint16_t DELTA> struct Unroller<DELTA, 31> {
  enum { SHIFT = (DELTA * 31) % 32 };

  static void Unpack(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
    out[31] = (*in) >> SHIFT ;
  }

  static void Pack(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
    *out |= (in[31] << SHIFT);
  }

  static void PackNoMask(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
    *out |= (in[31] << SHIFT);
  }

	// arcj: for iiu FIXME
  static void Unpack_iiu(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
		uint32_t shift = ((DELTA + 1) * 31) % 32;	
    out[31] = (*in) >> shift;
		//std::cout << "[arcj] Unpack_iiu2 " << DELTA << " / " << SHIFT << " / " << shift << " / " << in[31] << " / " << out[31] << std::endl;
  }

  static void Pack_iiu(const uint32_t *__restrict__ in, uint32_t *__restrict__ out, const uint32_t *__restrict__ exception, const size_t exceptcounter, const uint32_t count) {
		// case 3
		uint32_t shift = ((DELTA + 1) * 31) % 32;	
		//std::cout << "[arcj] Pack_iiu2 " << DELTA << " / " << SHIFT << " / " << shift << " / " << in[31] << std::endl;
		//bool exceptiontrue = false;
		//uint32_t pos = 0;
		//for(uint32_t z = 0; z < exceptcounter; ++z) {
		//	if(exception[z] == (31 + count)) {
		//		exceptiontrue = true;
		//		pos = z;
		//	}
		//}
		*out |= (in[31] << shift); 
		//if(exceptiontrue) {
		//	//std::cout << "[arcj] Pack_iiu2 exception-" << 31 + count << " / " << exception[pos] << " / " << pos << " >> " << pos + 128 << std::endl;
		//	*out |= ((pos + 128) << shift); // SHIFT == 8 
		//} else {
		//	//std::cout << "[arcj] Pack_iiu2 normal-" << 31 + count << " / " << in[31] << std::endl;
		//  *out |= (in[31] << shift); 
		//}
		//std::cout << "[arcj] OUT: " << *out << std::endl;
		//std::cout << "[arcj] out address3: " << out << std::endl;

    // *out |= (in[31] << SHIFT); // arcj: original
  }
};

}  // namespace


void __fastunpack0(const uint32_t *__restrict__, uint32_t *__restrict__ out) {
  for (uint32_t i = 0; i < 32; ++i)
    *(out++) = 0;
}
void __fastpack0(const uint32_t *__restrict__, uint32_t *__restrict__) {}
void __fastpackwithoutmask0(const uint32_t *__restrict__,
                            uint32_t *__restrict__) {}

void __fastunpack1(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<1>::Unpack(in, out);
}

void __fastunpack2(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<2>::Unpack(in, out);
}

void __fastunpack3(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<3>::Unpack(in, out);
}

void __fastunpack5(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<5>::Unpack(in, out);
}

void __fastunpack6(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<6>::Unpack(in, out);
}

void __fastunpack7(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<7>::Unpack(in, out);
}

void __fastunpack7_iiu(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<7>::Unpack_iiu(in, out);
}

void __fastunpack9(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  Unroller<9>::Unpack(in, out);
}

void __fastunpack10(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<10>::Unpack(in, out);
}

void __fastunpack11(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<11>::Unpack(in, out);
}

void __fastunpack12(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<12>::Unpack(in, out);
}

void __fastunpack13(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<13>::Unpack(in, out);
}

void __fastunpack14(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<14>::Unpack(in, out);
}

void __fastunpack15(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {

  Unroller<15>::Unpack(in, out);
}

void __fastunpack17(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<17>::Unpack(in, out);
}

void __fastunpack18(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<18>::Unpack(in, out);
}

void __fastunpack19(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<19>::Unpack(in, out);
}

void __fastunpack20(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<20>::Unpack(in, out);
}

void __fastunpack21(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<21>::Unpack(in, out);
}

void __fastunpack22(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<22>::Unpack(in, out);
}

void __fastunpack23(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<23>::Unpack(in, out);
}

void __fastunpack24(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<24>::Unpack(in, out);
}

void __fastunpack25(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<25>::Unpack(in, out);
}

void __fastunpack26(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<26>::Unpack(in, out);
}

void __fastunpack27(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<27>::Unpack(in, out);
}

void __fastunpack28(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<28>::Unpack(in, out);
}

void __fastunpack29(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<29>::Unpack(in, out);
}

void __fastunpack30(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<30>::Unpack(in, out);
}

void __fastunpack31(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  Unroller<31>::Unpack(in, out);
}

void __fastunpack32(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  for (int k = 0; k < 32; ++k)
    out[k] = in[k];
}
void __fastpack32(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  for (int k = 0; k < 32; ++k)
    out[k] = in[k];
}

void __fastpackwithoutmask32(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  for (int k = 0; k < 32; ++k)
    out[k] = in[k];
}

void __fastunpack4(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  for (uint32_t outer = 0; outer < 4; ++outer) {
    for (uint32_t inwordpointer = 0; inwordpointer < 32; inwordpointer += 4)
      *(out++) = ((*in) >> inwordpointer) % (1U << 4);
    ++in;
  }
}

void __fastunpack8(const uint32_t *__restrict__ in,
                   uint32_t *__restrict__ out) {
  for (uint32_t outer = 0; outer < 8; ++outer) {
    for (uint32_t inwordpointer = 0; inwordpointer < 32; inwordpointer += 8)
      *(out++) = ((*in) >> inwordpointer) % (1U << 8);
    ++in;
  }
}

void __fastunpack16(const uint32_t *__restrict__ in,
                    uint32_t *__restrict__ out) {
  for (uint32_t outer = 0; outer < 16; ++outer) {
    for (uint32_t inwordpointer = 0; inwordpointer < 32; inwordpointer += 16)
      *(out++) = ((*in) >> inwordpointer) % (1U << 16);
    ++in;
  }
}

void __fastpack1(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<1>::Pack(in, out);
}

void __fastpack2(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<2>::Pack(in, out);
}

void __fastpack3(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<3>::Pack(in, out);
}

void __fastpack5(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<5>::Pack(in, out);
}

void __fastpack6(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<6>::Pack(in, out);
}

void __fastpack7(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<7>::Pack(in, out);
}

void __fastpack7_iiu(const uint32_t *__restrict__ in, uint32_t *__restrict__ out, const uint32_t *__restrict__ exception, const size_t exceptcounter, const uint32_t count) {
  Unroller<7>::Pack_iiu(in, out, exception, exceptcounter, count);
}

void __fastpack9(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<9>::Pack(in, out);
}

void __fastpack10(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<10>::Pack(in, out);
}

void __fastpack11(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<11>::Pack(in, out);
}

void __fastpack12(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<12>::Pack(in, out);
}

void __fastpack13(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<13>::Pack(in, out);
}

void __fastpack14(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<14>::Pack(in, out);
}

void __fastpack15(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<15>::Pack(in, out);
}

void __fastpack17(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<17>::Pack(in, out);
}

void __fastpack18(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<18>::Pack(in, out);
}

void __fastpack19(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<19>::Pack(in, out);
}

void __fastpack20(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<20>::Pack(in, out);
}

void __fastpack21(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<21>::Pack(in, out);
}

void __fastpack22(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<22>::Pack(in, out);
}

void __fastpack23(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<23>::Pack(in, out);
}

void __fastpack24(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<24>::Pack(in, out);
}

void __fastpack25(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<25>::Pack(in, out);
}

void __fastpack26(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<26>::Pack(in, out);
}

void __fastpack27(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<27>::Pack(in, out);
}

void __fastpack28(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<28>::Pack(in, out);
}

void __fastpack29(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<29>::Pack(in, out);
}

void __fastpack30(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<30>::Pack(in, out);
}

void __fastpack31(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<31>::Pack(in, out);
}

void __fastpack4(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<4>::Pack(in, out);
}

void __fastpack8(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<8>::Pack(in, out);
}

void __fastpack16(const uint32_t *__restrict__ in, uint32_t *__restrict__ out) {
  Unroller<16>::Pack(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask1(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<1>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask2(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<2>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask3(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<3>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask5(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<5>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask6(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<6>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask7(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<7>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask9(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<9>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask10(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<10>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask11(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<11>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask12(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<12>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask13(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<13>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask14(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<14>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask15(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<15>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask17(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<17>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask18(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<18>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask19(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<19>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask20(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<20>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask21(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<21>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask22(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<22>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask23(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<23>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask24(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<24>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask25(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<25>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask26(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<26>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask27(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<27>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask28(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<28>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask29(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<29>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask30(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<30>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask31(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<31>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask4(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<4>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask8(const uint32_t *__restrict__ in,
                            uint32_t *__restrict__ out) {
  Unroller<8>::PackNoMask(in, out);
}

/*assumes that integers fit in the prescribed number of bits */
void __fastpackwithoutmask16(const uint32_t *__restrict__ in,
                             uint32_t *__restrict__ out) {
  Unroller<16>::PackNoMask(in, out);
}
