#include "util.h"

class fa_read {
public:
    fa_read(std::string id, size_t len, std::string seq_, std::string meta = "");
    void to_file(std::ostream& out);
    const size_t size() const;                                  // Returns read length
    const std::string& get_id() const;                          // Returns read id
    const std::string& get_seq() const;                         // Returns read
    const std::string& get_metadata() const;                    // Returns metadata iff it exists
    char operator[](size_t index) const;                        // Read only

protected:
    std::string read_id;
    std::string seq;
    std::string metadata;
    size_t length;
};