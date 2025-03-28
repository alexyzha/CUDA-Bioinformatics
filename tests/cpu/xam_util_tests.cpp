#include <gtest/gtest.h>
#include "../../src/cpu/headers/xam_util.h"
#include "../../src/cpu/headers/xam_parser.h"

TEST(XAM_UTIL, SAM_TO_FILE_SINGLE_NO_TAGS) {
    std::vector<sam_read*> reads = read_sam("infiles/sam_parser_read_no_tags_in.txt");
    std::string exp = "ID1\t0\tREF1\t1\tQ\tCIGAR\t=\t2\t4\tACTG\t????";
    EXPECT_EQ(reads.size(), 1) << RED << "GOT != 1 READ FROM FILE WITH 1 READ";
    std::ofstream out("outfiles/sam_to_file_no_tags_out.txt");
    EXPECT_TRUE(out.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    sam_read_to_file(out, *reads[0]);
    out.close();
    // Validate output
    std::ifstream in("outfiles/sam_to_file_no_tags_out.txt");
    EXPECT_TRUE(in.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    std::string line = "";
    std::getline(in, line);
    trim(line);
    EXPECT_EQ(line, exp) << RED << "MISMATCH BETWEEN EXP AND ACTUAL OUTPUT" << RESET << std::endl;
    EXPECT_FALSE(std::getline(in, line));
    in.close();
}

TEST(XAM_UTIL, SAM_TO_FILE_SINGLE_YES_TAGS) {
    std::vector<sam_read*> reads = read_sam("infiles/sam_parser_read_yes_tags_in.txt");
    std::string exp = "ID1\t0\tREF1\t1\tQ\tCIGAR\t=\t2\t4\tACTG\t????\tTAG1\tTAG2\tTAG3";
    EXPECT_EQ(reads.size(), 1) << RED << "GOT != 1 READ FROM FILE WITH 1 READ";
    std::ofstream out("outfiles/sam_to_file_yes_tags_out.txt");
    EXPECT_TRUE(out.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    sam_read_to_file(out, *reads[0]);
    out.close();
    // Validate output
    std::ifstream in("outfiles/sam_to_file_yes_tags_out.txt");
    EXPECT_TRUE(in.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    std::string line = "";
    std::getline(in, line);
    trim(line);
    EXPECT_EQ(line, exp) << RED << "MISMATCH BETWEEN EXP AND ACTUAL OUTPUT" << RESET << std::endl;
    EXPECT_FALSE(std::getline(in, line));
    in.close();
}

TEST(XAM_UTIL, SAM_TO_FILE_MULTI_ALL) {
    std::vector<sam_read*> reads = read_sam("infiles/sam_parser_all_allowed_in.txt");
    std::vector<std::string> exp = {
        "ID1\t1\tREF1\t10\tA\tAIGAR\t=\t100\t1000\tAAA\t???\tTAG1",
        "ID2\t2\tREF2\t20\tB\tBIGAR\t*\t200\t2000\tBBB\t?*!\tTAG2\tTAG2A",
        "ID3\t3\tREF3\t30\tC\tCIGAR\t!\t300\t3000\tCCC\t!!!\tTAG3\tTAG3A\tTAG3B"
    };
    EXPECT_EQ(reads.size(), 3) << RED << "GOT != 3 READ FROM FILE WITH 3 READS";
    std::ofstream out("outfiles/sam_to_file_multi_read_out.txt");
    EXPECT_TRUE(out.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    for(int i = 0; i < 3; ++i) {
        sam_read_to_file(out, *reads[i]);
    }
    out.close();
    // Validate output
    std::ifstream in("outfiles/sam_to_file_multi_read_out.txt");
    EXPECT_TRUE(in.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    std::string line = "";
    for(int i = 0; i < 3; ++i) {
        std::getline(in, line);
        trim(line);
        EXPECT_EQ(line, exp[i]) << RED << "MISMATCH BETWEEN EXP AND ACTUAL OUTPUT ON LINE [" << i << "]" << RESET << std::endl;
    }
    EXPECT_FALSE(std::getline(in, line));
    in.close();
}

TEST(XAM_UTIL, SAM_TO_FILE_MULTI_SORTED) {
    sam_container cont("infiles/sam_parser_all_allowed_in.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    std::vector<std::string> exp = {
        "ID3\t3\tREF3\t30\tC\tCIGAR\t!\t300\t3000\tCCC\t!!!\tTAG3\tTAG3A\tTAG3B",
        "ID2\t2\tREF2\t20\tB\tBIGAR\t*\t200\t2000\tBBB\t?*!\tTAG2\tTAG2A",
        "ID1\t1\tREF1\t10\tA\tAIGAR\t=\t100\t1000\tAAA\t???\tTAG1"
    };
    EXPECT_EQ(headers.size(), 3) << RED << "GOT != 3 HEADERS FROM FILE WITH 3 HEADERS" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 3) << RED << "GOT != 3 READS FROM FILE WITH 3 READS" << RESET << std::endl;
    cont.sort([&](sam_read* a, sam_read* b) {
        return a->flags > b->flags;
    });
    std::ofstream out("outfiles/sam_to_file_multi_sorted_out.txt");
    EXPECT_TRUE(out.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    for(int i = 0; i < 3; ++i) {
        sam_read_to_file(out, *reads[i]);
    }
    out.close();
    // Validate output
    std::ifstream in("outfiles/sam_to_file_multi_sorted_out.txt");
    EXPECT_TRUE(in.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    std::string line = "";
    for(int i = 0; i < 3; ++i) {
        std::getline(in, line);
        trim(line);
        EXPECT_EQ(line, exp[i]) << RED << "MISMATCH BETWEEN EXP AND ACTUAL OUTPUT ON LINE [" << i << "]" << RESET << std::endl;
    }
    EXPECT_FALSE(std::getline(in, line));
    in.close();
}

TEST(XAM_UTIL, SAM_CONTAINER_TO_FILE) {
    sam_container cont("infiles/sam_parser_all_allowed_in.txt");
    std::unordered_set<std::string> exp_header = {
        "@CA\tHEADER1",
        "@TZ\tHEADER2",
        "@CA\tHEADER3",
        "@TZ\tHEADER4",
        "@FG\tHEADER5"
    };
    std::vector<std::string> exp = {
        "ID1\t1\tREF1\t10\tA\tAIGAR\t=\t100\t1000\tAAA\t???\tTAG1",
        "ID2\t2\tREF2\t20\tB\tBIGAR\t*\t200\t2000\tBBB\t?*!\tTAG2\tTAG2A",
        "ID3\t3\tREF3\t30\tC\tCIGAR\t!\t300\t3000\tCCC\t!!!\tTAG3\tTAG3A\tTAG3B"
    };
    sam_to_file(cont, "outfiles/sam_container_out.txt");
    // Validate output
    std::ifstream in("outfiles/sam_container_out.txt");
    EXPECT_TRUE(in.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    std::string line = "";
    for(int i = 0; i < 5; ++i) {
        std::getline(in, line);
        trim(line);
        EXPECT_TRUE(exp_header.count(line)) << RED << "UNEXPECTED HEADER IN OUTFILE" << RESET << std::endl;
        exp_header.erase(line);
    }
    EXPECT_TRUE(exp_header.empty()) << RED << "MISSING HEADER(S) IN OUTPUT" << RESET << std::endl;
    for(int i = 0; i < 3; ++i) {
        std::getline(in, line);
        trim(line);
        EXPECT_EQ(line, exp[i]) << RED << "MISMATCH BETWEEN EXP AND ACTUAL OUTPUT ON LINE [" << i + 5 << "]" << RESET << std::endl;
    }
    EXPECT_FALSE(std::getline(in, line));
    in.close();
}

TEST(XAM_UTIL, SAM_CONTAINER_TO_FILE_SORTED_READS) {
    sam_container cont("infiles/sam_parser_all_allowed_in.txt");
    std::unordered_set<std::string> exp_header = {
        "@CA\tHEADER1",
        "@TZ\tHEADER2",
        "@CA\tHEADER3",
        "@TZ\tHEADER4",
        "@FG\tHEADER5"
    };
    std::vector<std::string> exp = {
        "ID3\t3\tREF3\t30\tC\tCIGAR\t!\t300\t3000\tCCC\t!!!\tTAG3\tTAG3A\tTAG3B",
        "ID2\t2\tREF2\t20\tB\tBIGAR\t*\t200\t2000\tBBB\t?*!\tTAG2\tTAG2A",
        "ID1\t1\tREF1\t10\tA\tAIGAR\t=\t100\t1000\tAAA\t???\tTAG1"
    };
    cont.sort([&](sam_read* a, sam_read* b){
        return a->qname > b->qname;
    });
    sam_to_file(cont, "outfiles/sam_container_sorted_out.txt");
    // Validate output
    std::ifstream in("outfiles/sam_container_sorted_out.txt");
    EXPECT_TRUE(in.is_open()) << RED << "UNABLE TO OPEN FILE" << RESET << std::endl;
    std::string line = "";
    for(int i = 0; i < 5; ++i) {
        std::getline(in, line);
        trim(line);
        EXPECT_TRUE(exp_header.count(line)) << RED << "UNEXPECTED HEADER IN OUTFILE" << RESET << std::endl;
        exp_header.erase(line);
    }
    EXPECT_TRUE(exp_header.empty()) << RED << "MISSING HEADER(S) IN OUTPUT" << RESET << std::endl;
    for(int i = 0; i < 3; ++i) {
        std::getline(in, line);
        trim(line);
        EXPECT_EQ(line, exp[i]) << RED << "MISMATCH BETWEEN EXP AND ACTUAL OUTPUT ON LINE [" << i + 5 << "]" << RESET << std::endl;
    }
    EXPECT_FALSE(std::getline(in, line));
    in.close();
}

TEST(XAM_UTIL, CIGAR_COMPLETE_MATCH) {
    std::string* ref = new std::string("AAAA");
    std::string* read = new std::string("AAAA");
    std::string exp = "4M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 4M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_MISMATCH) {
    std::string* ref = new std::string("AAAGA");
    std::string* read = new std::string("AAAAA");
    std::string exp = "3M1X1M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 3M1X1M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_REF_GAP) {
    std::string* ref = new std::string("AAA-A");
    std::string* read = new std::string("AAAAA");
    std::string exp = "3M1I1M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 3M1I1M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_READ_GAP) {
    std::string* ref = new std::string("AAAAA");
    std::string* read = new std::string("AA-AA");
    std::string exp = "2M1D2M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 2M1D2M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_MULTI_MISMATCH) {
    std::string* ref = new std::string("ACTGA");
    std::string* read = new std::string("AAAAA");
    std::string exp = "1M3X1M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 1M3X1M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_MULTI_REF_GAP) {
    std::string* ref = new std::string("AA--A");
    std::string* read = new std::string("AAAAA");
    std::string exp = "2M2I1M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 2M2I1M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_MULTI_READ_GAP) {
    std::string* ref = new std::string("AAAAA");
    std::string* read = new std::string("A--AA");
    std::string exp = "1M2D2M";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 1M2D2M, NOT [" << cigar << "]" << RESET << std::endl; 
}

TEST(XAM_UTIL, CIGAR_ALL_ALLOWED) {
    std::string* ref = new std::string("ACG---CGT--ACG-TAGC");
    std::string* read = new std::string("GCGACTT--ACACGATA--");
    std::string exp = "1X2M3I1X2D2I3M1I2M2D";
    alignment align(0, 0, 0, ref, read);
    std::string cigar = make_cigar(align);
    EXPECT_EQ(cigar, exp) << RED << "CIGAR SHOULD = 1X2M3I1X2D2I3M1I2M2D, NOT [" << cigar << "]" << RESET << std::endl;
}

TEST(XAM_UTIL, MAP_READS_TO_REF_MATCH_SINGLE) {
    std::string ref = "AAAA";
    std::string ref_id = "REFID";
    std::string read = "AAAA";
    std::vector<fq_read*> fqreads = {};
    fqreads.push_back(new fq_read("READID", 4, read, "!!!!", ""));
    std::vector<sam_read*> samreads = map_reads_to_ref(ref, ref_id, fqreads, 4);
    // Validate samreads
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*samreads[0]);
    EXPECT_TRUE(tags.empty()) << RED << "NO TAGS SHOULD BE GENERATED" << RESET << std::endl;
    EXPECT_EQ(qname, "READID") << RED << "DIFFERENT QNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rname, ref_id) << RED << "DIFFERENT RNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(cigar, "4M") << RED << "DIFFERENT CIGAR THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rnext, "*") << RED << "DIFFERENT RNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(seq, read) << RED << "DIFFERENT SEQ THAN EXP" << RESET << std::endl;
    EXPECT_EQ(qual, "!!!!") << RED << "DIFFERENT QUAL THAN EXP" << RESET << std::endl;
    EXPECT_EQ(flags, 0) << RED << "DIFFERENT FLAGS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(pos, 0) << RED << "DIFFERENT POS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(posnext, 0) << RED << "DIFFERENT POSNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(tlen, 0) << RED << "DIFFERENT TLEN THAN EXP" << RESET << std::endl;
    EXPECT_EQ(mapq, static_cast<char>(255)) << RED << "DIFFERENT MAPQ THAN EXP" << RESET << std::endl;
}

TEST(XAM_UTIL, MAP_READS_TO_REF_MIDDLE_SINGLE) {
    std::string ref = "AAGCTAAA";
    std::string ref_id = "REFID";
    std::string read = "GCT";
    std::vector<fq_read*> fqreads = {};
    fqreads.push_back(new fq_read("READID", 3, read, "!!!", ""));
    std::vector<sam_read*> samreads = map_reads_to_ref(ref, ref_id, fqreads, 3);
    // Validate samreads
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*samreads[0]);
    EXPECT_TRUE(tags.empty()) << RED << "NO TAGS SHOULD BE GENERATED" << RESET << std::endl;
    EXPECT_EQ(qname, "READID") << RED << "DIFFERENT QNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rname, ref_id) << RED << "DIFFERENT RNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(cigar, "3M") << RED << "DIFFERENT CIGAR THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rnext, "*") << RED << "DIFFERENT RNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(seq, read) << RED << "DIFFERENT SEQ THAN EXP" << RESET << std::endl;
    EXPECT_EQ(qual, "!!!") << RED << "DIFFERENT QUAL THAN EXP" << RESET << std::endl;
    EXPECT_EQ(flags, 0) << RED << "DIFFERENT FLAGS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(pos, 2) << RED << "DIFFERENT POS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(posnext, 0) << RED << "DIFFERENT POSNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(tlen, 0) << RED << "DIFFERENT TLEN THAN EXP" << RESET << std::endl;
    EXPECT_EQ(mapq, static_cast<char>(255)) << RED << "DIFFERENT MAPQ THAN EXP" << RESET << std::endl;
}

TEST(XAM_UTIL, MAP_READS_TO_REF_MIDDLE_MISMATCH_SINGLE) {
    std::string ref = "AAGCTAAA";
    std::string ref_id = "REFID";
    std::string read = "GAT";
    std::vector<fq_read*> fqreads = {};
    fqreads.push_back(new fq_read("READID", 3, read, "!!!", ""));
    std::vector<sam_read*> samreads = map_reads_to_ref(ref, ref_id, fqreads, 3);
    // Validate samreads
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*samreads[0]);
    EXPECT_TRUE(tags.empty()) << RED << "NO TAGS SHOULD BE GENERATED" << RESET << std::endl;
    EXPECT_EQ(qname, "READID") << RED << "DIFFERENT QNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rname, ref_id) << RED << "DIFFERENT RNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(cigar, "*") << RED << "DIFFERENT CIGAR THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rnext, "*") << RED << "DIFFERENT RNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(seq, read) << RED << "DIFFERENT SEQ THAN EXP" << RESET << std::endl;
    EXPECT_EQ(qual, "!!!") << RED << "DIFFERENT QUAL THAN EXP" << RESET << std::endl;
    EXPECT_EQ(flags, 0x4) << RED << "DIFFERENT FLAGS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(pos, 0) << RED << "DIFFERENT POS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(posnext, 0) << RED << "DIFFERENT POSNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(tlen, 0) << RED << "DIFFERENT TLEN THAN EXP" << RESET << std::endl;
    EXPECT_EQ(mapq, static_cast<char>(255)) << RED << "DIFFERENT MAPQ THAN EXP" << RESET << std::endl;
}

TEST(XAM_UTIL, MAP_READS_TO_REF_MIDDLE_MISMATCH_SINGLE_K_ALLOW) {
    std::string ref = "AAGCTAAA";
    std::string ref_id = "REFID";
    std::string read = "GCTCA";
    std::vector<fq_read*> fqreads = {};
    fqreads.push_back(new fq_read("READID", 5, read, "!!!!!", ""));
    std::vector<sam_read*> samreads = map_reads_to_ref(ref, ref_id, fqreads, 3);
    // Validate samreads
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*samreads[0]);
    EXPECT_TRUE(tags.empty()) << RED << "NO TAGS SHOULD BE GENERATED" << RESET << std::endl;
    EXPECT_EQ(qname, "READID") << RED << "DIFFERENT QNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rname, ref_id) << RED << "DIFFERENT RNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(cigar, "3M1X1M") << RED << "DIFFERENT CIGAR THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rnext, "*") << RED << "DIFFERENT RNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(seq, read) << RED << "DIFFERENT SEQ THAN EXP" << RESET << std::endl;
    EXPECT_EQ(qual, "!!!!!") << RED << "DIFFERENT QUAL THAN EXP" << RESET << std::endl;
    EXPECT_EQ(flags, 0) << RED << "DIFFERENT FLAGS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(pos, 2) << RED << "DIFFERENT POS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(posnext, 0) << RED << "DIFFERENT POSNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(tlen, 0) << RED << "DIFFERENT TLEN THAN EXP" << RESET << std::endl;
    EXPECT_EQ(mapq, static_cast<char>(255)) << RED << "DIFFERENT MAPQ THAN EXP" << RESET << std::endl;
}

TEST(XAM_UTIL, MAP_READS_TO_REF_MIDDLE_GAP_SINGLE_K_ALLOW) {
    std::string ref = "AAGCTACCC";
    std::string ref_id = "REFID";
    std::string read = "GCTCCC";
    std::vector<fq_read*> fqreads = {};
    fqreads.push_back(new fq_read("READID", 6, read, "!!!!!!", ""));
    std::vector<sam_read*> samreads = map_reads_to_ref(ref, ref_id, fqreads, 3);
    // Validate samreads
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*samreads[0]);
    EXPECT_TRUE(tags.empty()) << RED << "NO TAGS SHOULD BE GENERATED" << RESET << std::endl;
    EXPECT_EQ(qname, "READID") << RED << "DIFFERENT QNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rname, ref_id) << RED << "DIFFERENT RNAME THAN EXP" << RESET << std::endl;
    EXPECT_EQ(cigar, "3M1D3M") << RED << "DIFFERENT CIGAR THAN EXP" << RESET << std::endl;
    EXPECT_EQ(rnext, "*") << RED << "DIFFERENT RNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(seq, read) << RED << "DIFFERENT SEQ THAN EXP" << RESET << std::endl;
    EXPECT_EQ(qual, "!!!!!!") << RED << "DIFFERENT QUAL THAN EXP" << RESET << std::endl;
    EXPECT_EQ(flags, 0) << RED << "DIFFERENT FLAGS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(pos, 2) << RED << "DIFFERENT POS THAN EXP" << RESET << std::endl;
    EXPECT_EQ(posnext, 0) << RED << "DIFFERENT POSNEXT THAN EXP" << RESET << std::endl;
    EXPECT_EQ(tlen, 0) << RED << "DIFFERENT TLEN THAN EXP" << RESET << std::endl;
    EXPECT_EQ(mapq, static_cast<char>(255)) << RED << "DIFFERENT MAPQ THAN EXP" << RESET << std::endl;
}

TEST(XAM_UTIL, MAP_READS_TO_REF_MATCH_MULTI) {
    std::string ref = "CCGTACACTG";
    std::string ref_id = "REFID";
    std::vector<std::string> reads = {"GTA", "ACA", "CTG"};
    std::string exp_cigar = "3M";
    std::vector<int> exp_pos = {2, 4, 7};
    std::vector<fq_read*> fqreads = {};
    for(int i = 0; i < 3; ++i) {
        fqreads.push_back(new fq_read("READID" + std::to_string(i), 3, reads[i], "!!!", ""));
    }
    std::vector<sam_read*> samreads = map_reads_to_ref(ref, ref_id, fqreads, 3);
    // Validate samreads
    for(int i = 0; i < 3; ++i) {
        auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*samreads[i]);
        EXPECT_TRUE(tags.empty()) << RED << "NO TAGS SHOULD BE GENERATED" << RESET << std::endl;
        EXPECT_EQ(qname, "READID" + std::to_string(i)) << RED << "DIFFERENT QNAME THAN EXP" << RESET << std::endl;
        EXPECT_EQ(rname, ref_id) << RED << "DIFFERENT RNAME THAN EXP" << RESET << std::endl;
        EXPECT_EQ(cigar, exp_cigar) << RED << "DIFFERENT CIGAR THAN EXP" << RESET << std::endl;
        EXPECT_EQ(rnext, "*") << RED << "DIFFERENT RNEXT THAN EXP" << RESET << std::endl;
        EXPECT_EQ(seq, reads[i]) << RED << "DIFFERENT SEQ THAN EXP" << RESET << std::endl;
        EXPECT_EQ(qual, "!!!") << RED << "DIFFERENT QUAL THAN EXP" << RESET << std::endl;
        EXPECT_EQ(flags, 0) << RED << "DIFFERENT FLAGS THAN EXP" << RESET << std::endl;
        EXPECT_EQ(pos, exp_pos[i]) << RED << "DIFFERENT POS THAN EXP" << RESET << std::endl;
        EXPECT_EQ(posnext, 0) << RED << "DIFFERENT POSNEXT THAN EXP" << RESET << std::endl;
        EXPECT_EQ(tlen, 0) << RED << "DIFFERENT TLEN THAN EXP" << RESET << std::endl;
        EXPECT_EQ(mapq, static_cast<char>(255)) << RED << "DIFFERENT MAPQ THAN EXP" << RESET << std::endl;
    }
}

TEST(XAM_UTIL, SAM_TO_VCF_EMPTY) {
    // Also technically "SAM_TO_VCF_HEADER_ONLY"
    std::vector<sam_read*> samreads;
    std::string ref_seq = "REFSEQ";
    std::string line = "";
    std::vector<std::string> exp_headers = {
        "##fileformat=VCFv4.2",
        "##source=cubio",
        "##reference=REFID",
        "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    };
    sam_to_vcf("outfiles/sam_to_vcf_empty_out.txt", ref_seq, "REFID", samreads, 7);
    // Validate output
    std::ifstream file("outfiles/sam_to_vcf_empty_out.txt");
    EXPECT_TRUE(file.is_open()) << RED << "INVALID FILE OR FILE DOESN'T EXIST" << RESET << std::endl;
    for(int i = 0; i < 4; ++i) {
        std::getline(file, line);
        trim(line);
        EXPECT_EQ(line, exp_headers[i]) << RED << "WRONG HEADER ON LINE [" << i << "]" << RESET << std::endl;
    }
    EXPECT_FALSE(std::getline(file, line)) << RED << "EXTRA LINES IN FILE" << RESET << std::endl;
    file.close();
}

TEST(XAM_UTIL, SAM_TO_VCF_ALL) {
    std::vector<sam_read*> samreads = read_sam("infiles/sam_to_vcf_all_in.txt");
    std::string ref_seq = "ACTT";
    std::string line = "";
    std::vector<std::string> exp_lines = {
        "##fileformat=VCFv4.2",
        "##source=cubio",
        "##reference=REFID",
        "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "7\t4\t.\tT\tG\t255\tPASS\tDP=1",
        "7\t2\t.\tC\tG\t255\tPASS\tDP=1",
        "7\t1\t.\tA\tC\t255\tPASS\tDP=1"
    };
    sam_to_vcf("outfiles/sam_to_vcf_all_out.txt", ref_seq, "REFID", samreads, 7);
    // Validate output
    std::ifstream file("outfiles/sam_to_vcf_all_out.txt");
    EXPECT_TRUE(file.is_open()) << RED << "INVALID FILE OR FILE DOESN'T EXIST" << RESET << std::endl;
    for(int i = 0; i < 7; ++i) {
        std::getline(file, line);
        trim(line);
        EXPECT_EQ(line, exp_lines[i]) << RED << "WRONG VALUE ON LINE [" << i << "]" << RESET << std::endl;
    }
    EXPECT_FALSE(std::getline(file, line)) << RED << "EXTRA LINES IN FILE" << RESET << std::endl;
    file.close();
}