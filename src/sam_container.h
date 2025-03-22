#include "sam_read.h"
#include "cpu/xam_parser.h"

class sam_container {
public:
    sam_container(std::string file_path = "");
    const std::unordered_map<std::string, std::vector<std::string>>& get_headers() const;
    const std::vector<sam_read*>& get_reads() const;
    void set_headers(std::string file_path);
    void set_reads(std::string file_path);

    template<typename T>
    void sort(T comp) {
        std::sort(reads.begin(), reads.end(), comp);
    }

private:
    std::vector<sam_read*> reads;
    std::unordered_map<std::string, std::vector<std::string>> headers;
    
};