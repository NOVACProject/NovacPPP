#pragma once

#include <cstdint>

namespace SpectrumIO
{

	#define hdr_version 5

	struct MKZYhdr
	{
		char ident[4];                  // "MKZY"
		unsigned short hdrsize;         // this is the size in bytes of the header
		unsigned short hdrversion;      // version of the header
		unsigned short size;            // the number of bytes with compressed data
		unsigned short checksum;        // checksum for the uncompressed data
		char name[12];                  // the name of this specific measurement
		char instrumentname[16];        // the name of the instrument
		unsigned short startc;          // the startchannel for the first data-point
		unsigned short pixels;          // number of pixels saved in the data-field
		short viewangle;                // the viewing angle of the instrument
		unsigned short scans;           // total number of scans added
		short exptime;                  // exposure time, negative if set automatic
		unsigned char channel;          // channel of the spectrometer, typically 0
		unsigned char flag;             // for further use, currently contains the
										// status of the solenoid(s) in bit 0 and 1
		std::uint32_t date;             // date. Type changed from 'unsigned long' to uint32_t on 2018-01-13 by Mattias
		std::uint32_t starttime;        // time when the scanning was started. Type changed from 'unsigned long' to uint32_t on 2018-01-13 by Mattias
		std::uint32_t stoptime;         // time when the scanning was finished. Type changed from 'unsigned long' to uint32_t on 2018-01-13 by Mattias
		double lat;                     // GPS latitude in degrees
		double lon;                     // GPS longitude in degrees
		short altitude;                 // new in version 2
		char measureidx;                // new in version 2, nr between 0 and measurecnt-1
		char measurecnt;                 //new in version 2
										// number of MEAS= lines in cfg.txt
		short viewangle2;                //new in version 3, direction of 2nd motor
		short compassdir;                //new in version 3, given in cfg.txt
		short tiltX;                     //new in version 3, given in cfg.txt
		short tiltY;                     //new in version 3, given in cfg.txt
		float temperature;               //new in version 3, given in cfg.txt
		char coneangle;                  //new in version 4, given in cfg.txt
		unsigned short ADC[8];           //new in version 5
	};

	#define headsiz 12

	/** <b>MKPack</b> is a class for reading and writing spectra in 
		Manne Kihlman's Pak-format. */
	class MKPack
	{
	public:
		MKPack();
		~MKPack();

		/** Compress the given input-buffer to the given output-buffer 
			@param in - the uncompressed spectral data
			@param ut - will on return contain the compressed spectrum
			@param size - the number of data-points in the uncompressed spectrum
			@return - the number of data-points in the compressed spectrum */
		unsigned short mk_compress(long *in,unsigned char *ut,unsigned short size);

		/** Uncompress the given input-buffer to the given output-buffer
				@param inpek - the input-buffer
				@param kvar - the number of data-points to be written to the output-buffer
				@param ut - will in return contain the uncompressed spectrum
				@return -  the number of data-points in the uncompressed spectrum ???? */
		long UnPack(unsigned char *inpek, long kvar, long *ut );

	private:
		void SetBit(unsigned char *pek, long bit);
		void ClearBit(unsigned char *pek,long bit);

		void WriteBits(short a,short curr,long *inpek,unsigned char *utpek,long bitnr);
		short BitsPrec(long i);

		void PackSeg(unsigned char *utpek, long *kvar);

		long m_bitnr;
		long *m_strt;

	};
}