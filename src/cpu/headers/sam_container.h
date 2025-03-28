#include "sam_read.h"
#include "xam_parser.h"

class sam_container {
public:
    /*
     *  `sam_container` constructor. Contains `sam_read`.
     *  @param file_path `string` reads sam reads and headers from file, default = ""
     */
    sam_container(std::string file_path = "");

    /*
     *  Returns stored headers.
     *  @return `const unordered_map<string, vector<string>>&`
     */
    const std::unordered_map<std::string, std::vector<std::string>>& get_headers() const;

    /*
     *  Returns stored sam reads.
     *  @return `const vector<sam_read*>&`
     */
    const std::vector<sam_read*>& get_reads() const;

    /*
     *  Sets stored headers from file.
     *  @param file_path `string` reads headers from file
     *  @return `void`
     */
    void set_headers(std::string file_path);

    /*
     *  Sets stored sam reads from file.
     *  @param file_path `string` reads sam reads from file
     *  @return `void`
     */
    void set_reads(std::string file_path);

    /*
     *  Allows for a custom sorting scheme for stored sam reads.
     *  @param comp `*` comparator for sorting sam reads
     *  @param comp decl ex: [&](sam_read* a, sam_read* b) {...}
     *  @return `void`
     */
    template<typename T>
    void sort(T comp) {
        std::sort(reads.begin(), reads.end(), comp);
    }

private:
    std::vector<sam_read*> reads;
    std::unordered_map<std::string, std::vector<std::string>> headers;
    
};