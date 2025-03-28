#include "fa_read.h"

#ifndef C_FQ_READ
#define C_FQ_READ

class fq_read : public fa_read {
public:
    /*
     *  `fq_read` constructor. Inherits from `fa_read`.
     *  @param id `string` sequence id
     *  @param len `size_t` length of sequence
     *  @param seq_ `string` sequence
     *  @param qual_ `string` quality of sequence
     *  @param meta `string` sequence metadata, default = ""
     */
    fq_read(std::string id, size_t len, std::string seq_, std::string qual_, std::string meta = "");
    
    /*
     *  Streams fasta read to a file.
     *  @param out `ofstream&` file out stream handle to print contents to 
     *  @return `void`
     */
    void to_file(std::ofstream& out);

    /*
     *  Returns reference to quality of fastq read.
     *  @return `string&`
     */
    const std::string& get_quality() const;

    /*
     *  Returns 16-bit mask with both sequence[index] and quality[index] of fastq read.
     *  sequence[index] is in bits `0-7`, quality[index] is in bits `8-15`.
     *  @param index `size_t` desired index of sequence
     *  @return `uint16_t` EXPECT(0 <= index <= sequence length)
     */
    uint16_t operator[](size_t index) const;
    
private:
    std::string qual;
    
};

#endif