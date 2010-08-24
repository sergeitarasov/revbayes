#ifndef RandomNumberFactory_H
#define RandomNumberFactory_H

#include <set>
#include <vector>

#define GLOBAL_RNG RandomNumberFactory::randomNumberFactoryInstance().getGlobalRandomNumberGenerator()
#define NEW_RNG RandomNumberFactory::randomNumberFactoryInstance().getRandomNumberGenerator()

class RandomNumberGenerator;
class RandomNumberFactory {


public:
    static RandomNumberFactory& randomNumberFactoryInstance(void)
        {
        static RandomNumberFactory singleRandomNumberFactory;
        return singleRandomNumberFactory;
        }
    void deleteRandomNumberGenerator(RandomNumberGenerator* r);
    RandomNumberGenerator* getRandomNumberGenerator(void);
    RandomNumberGenerator* getRandomNumberGenerator(std::vector<unsigned int> s);
    RandomNumberGenerator* getGlobalRandomNumberGenerator(void) { return seedGenerator; }
    void setFactorySeed(unsigned int &seed1, unsigned int &seed2) { I1 = seed1; I2 = seed2; }

private:
    RandomNumberFactory(void);
    RandomNumberFactory(const RandomNumberFactory&);
    RandomNumberFactory&operator=(const RandomNumberFactory&);
    ~RandomNumberFactory(void);
    void setSeed(void);
    RandomNumberGenerator *seedGenerator;
    std::set<RandomNumberGenerator*> allocatedRandomNumbers;
    unsigned int I1, I2;
};

#endif


