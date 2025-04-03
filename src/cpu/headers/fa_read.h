#include "util.h"

#ifndef C_FA_READ
#define C_FA_READ

/**
 *  @brief Fasta read container.
 */
class fa_read {
public:
    /**
     *  @brief `fa_read` constructor.
     *  @param id `string` sequence id
     *  @param len `size_t` length of sequence
     *  @param seq_ `string` sequence
     *  @param meta `string` sequence metadata, default = ""
     *  @note Copies string inputs.
     */
    fa_read(std::string id, size_t len, std::string seq_, std::string meta = "");

    /**
     *  @brief Streams fasta read to a file.
     *  @param out `ofstream&` file out stream handle to print contents to 
     *  @return `void`
     */
    void to_file(std::ofstream& out);

    /**
     *  @brief Returns length of fasta read sequence.
     *  @return `size_t`
     */
    const size_t size() const;

    /**
     *  @brief Returns reference to ID of fasta read.
     *  @return `string&`
     */
    const std::string& get_id() const;

    /**
     *  @brief Returns reference to fasta read sequence.
     *  @return `string&`
     */
    const std::string& get_seq() const;

    /**
     *  @brief Returns reference to fasta read metadata.
     *  @return `string&`
     */
    const std::string& get_metadata() const;
    
    /**
     *  @brief Returns fasta read sequence[index].
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