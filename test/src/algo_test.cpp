#include "test_common.hpp"

class Algo : public ::testing::Test
{
public:
    std::string path = "genome_samples/s2286.txt";
    std::vector<uint32_t>freq1;
    std::vector<uint32_t>freq2;
    
};

TEST_F(Algo, test_naive_scalar)
{
    uint32_t lenght = rand()%1020 + 4;
    naive_scalar(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_hash3_scalar)
{
    uint32_t lenght = rand()%1020 + 4;
    
    hash3_scalar(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_rabin_karp_rolling_hash_scalar)
{
    uint32_t lenght = rand()%100 + 10;
    rabin_karp_rolling_hash_scalar(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_rabin_karp_SWAR_scalar)
{
    uint32_t lenght = rand()%100 + 10;
    rabin_karp_SWAR_scalar(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_naive_vector)
{
    uint32_t lenght = rand()%1020 + 4;
    naive_vector(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_hash3_vector)
{
    uint32_t lenght = rand()%1020 + 4;
    
    hash3_vector(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_rabin_karp_rolling_hash_vector)
{
    uint32_t lenght = rand()%100 + 10;
    rabin_karp_rolling_hash_vector(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}

TEST_F(Algo, test_rabin_karp_SWAR_vector)
{
    uint32_t lenght = rand()%100 + 10;
    rabin_karp_SWAR_vector(freq1, path, lenght,false);
    base_naive(freq2, path, lenght);
    EXPECT_EQ(freq1,freq2);
}