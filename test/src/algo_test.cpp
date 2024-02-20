#include "test_common.hpp"

class Algo : public ::testing::Test
{
public:
    std::string path = "genome_samples/s2286.txt";
    std::vector<uint32_t>freq1;
    std::vector<uint32_t>freq2;
    
};

TEST_F(Algo, test_naive)
{
    uint32_t lenght = rand()%1020 + 4;
    naive(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_hash3)
{
    uint32_t lenght = rand()%1020 + 4;
    
    hash3(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_rabin_karp)
{
    uint32_t lenght = rand()%1020 + 4;
    rabin_karp(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}