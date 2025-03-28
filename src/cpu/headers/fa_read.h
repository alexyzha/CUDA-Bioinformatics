#include "util.h"

#ifndef C_FA_READ
#define C_FA_READ

class fa_read {
public:
    /*
     *  `fa_read` constructor.
     *  @param id `string` sequence id
     *  @param len `size_t` length of sequence
     *  @param seq_ `string` sequence
     *  @param meta `string` sequence metadata, default = ""
     */
    fa_read(std::string id, size_t len, std::string seq_, std::string meta = "");

    /*
     *  Streams fasta read to a file.
     *  @param out `ofstream&` file out stream handle to print contents to 
     *  @return `void`
     */
    void to_file(std::ofstream& out);

    /*
     *  Returns length of fasta read sequence.
     *  @return `size_t`
     */
    const size_t size() const;

    /*
     *  Returns reference to ID of fasta read.
     *  @return `string&`
     */
    const std::string& get_id() const;

    /*
     *  Returns reference to fasta read sequence.
     *  @return `string&`
     */
    const std::string& get_seq() const;

    /*
     *  Returns reference to fasta read metadata.
     *  @return `string&`
     */
    const std::string& get_metadata() const;
    
    /*
     *  Returns fasta read sequence[index].
     *  @param index `size_t` desired index of sequence
     *  @return `char` EXPECT(0 <= index <= sequence)
     */
    char operator[](size_t index) const;

protected:
    std::string read_id;
    std::string seq;
    std::string metadata;
    size_t length;
    
};

#endif