// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_MUTATIONPROCESS_H
#define BPP_PHYL_SIMULATION_MUTATIONPROCESS_H

#include <Bpp/Numeric/VectorTools.h>

#include "../Mapping/SubstitutionRegister.h"
#include "../Model/SubstitutionModel.h"

namespace bpp
{
/**
 * @brief This class is used by MutationProcess to store detailed results of simulations.
 *
 * @author Julien Dutheil
 */
class MutationPath
{
private:
  std::shared_ptr<const Alphabet> alphabet_;

  /**
   * @brief The states taken, without initial state.
   */
  std::vector<size_t> states_;

  /**
   * @brief Times between states.
   * The first element in array is the time between the initial state and the first state in states_.
   */
  std::vector<double> times_;

  /**
   * @brief The initial state.
   */
  size_t initialState_;

  /**
   * @brief Total time of evolution.
   * Typically, this is a branch length.
   */
  double totalTime_;

public:
  /**
   * @brief Builds a new MutationPath object with initial state 'initialState' and total time 'time'.
   *
   * @param alphabet     The alphabet associated to the states in this path.
   * @param initialState The initial state.
   * @param time         The total time of evolution.
   */
  MutationPath(
      std::shared_ptr<const Alphabet> alphabet,
      size_t initialState,
      double time) :
    alphabet_(alphabet), states_(), times_(), initialState_(initialState), totalTime_(time) {}

  MutationPath(const MutationPath& path) :
    alphabet_(path.alphabet_), states_(path.states_), times_(path.times_), initialState_(path.initialState_), totalTime_(path.totalTime_) {}

  MutationPath& operator=(const MutationPath& path)
  {
    alphabet_     = path.alphabet_;
    states_       = path.states_;
    times_        = path.times_;
    initialState_ = path.initialState_;
    totalTime_    = path.totalTime_;
    return *this;
  }

  virtual ~MutationPath() {}

public:
  /**
   * @return A pointer toward the alphabet associated to this path.
   */
  std::shared_ptr<const Alphabet> getAlphabet() const { return alphabet_; }

  /**
   * @return A reference toward the alphabet associated to this path.
   */
  const Alphabet& alphabet() const { return *alphabet_; }

  /**
   * @brief Add a new mutation event.
   *
   * @param state The new state after mutation event.
   * @param time  The time between this mutation and previous mutation (or initial state).
   */
  void addEvent(size_t state, double time)
  {
    states_.push_back(state);
    times_.push_back(time);
  }

  /*
   * @brief Remove all mutations
   *
   */
  void clear()
  {
    states_.clear();
    times_.clear();
  }

  /**
   * @brief Retrieve the initial state.
   *
   * @return The initial state of this path.
   */
  size_t getInitialState() const { return initialState_; }

  /**
   * @brief Retrieve the total time of evolution.
   *
   * @return The total time of evolution.
   */
  double getTotalTime() const { return totalTime_; }

  /**
   * @brief Retrieve the number of substitution events.
   *
   * @return The number of substitution events, i.e. the number of states (without initial state).
   */
  size_t getNumberOfEvents() const { return states_.size(); }

  /**
   * @brief Retrieve the number of substitution events per type of substitution.
   *
   * @param counts A matrix with the same size as the alphabet. The substitution counts will be incremented according to the mutation path, which allows to efficiently sum various mutation paths with a look.
   */
  template<class Scalar>
  void getEventCounts(Matrix<Scalar>& counts) const
  {
    if (counts.getNumberOfRows()    != alphabet_->getSize()
        || counts.getNumberOfColumns() != alphabet_->getSize())
      throw Exception("MutationPath::getEventCounts. Incorrect input matrix, does not match alphabet size.");
    size_t currentState = initialState_;
    for (size_t i = 0; i < states_.size(); ++i)
    {
      size_t newState = states_[i];
      counts(currentState, newState)++;
      currentState = newState;
    }
  }

  /**
   * @brief Retrieve the number of substitution events per type of substitution, defined by a SubstitutionRegister object.
   *
   * @param counts A vector with the appropriate size, as defined by SubstitutionRegister::getNumberOfSubstitutionTypes(). The substitution counts will be incremented according to the mutation path, which allows to efficiently sum various mutation paths with a look.
   * @param reg The substitution register to use to categorize substitutions.
   */
  template<class Scalar>
  void getEventCounts(std::vector<Scalar>& counts, const SubstitutionRegisterInterface& reg) const
  {
    if (counts.size() != reg.getNumberOfSubstitutionTypes())
      throw Exception("MutationPath::getEventCounts. Incorrect input vector, does not match alphabet size.");
    size_t currentState = initialState_;
    for (size_t i = 0; i < states_.size(); ++i)
    {
      size_t newState = states_[i];
      size_t type = reg.getType(currentState, newState);
      if (type > 0) counts[type - 1]++;
      currentState = newState;
    }
  }

  /**
   * @brief Retrieve the final state of this path.
   *
   * @return The initial state if no mutation occurred, otherwise sends the state after last mutation event.
   */
  size_t getFinalState() const
  {
    if (states_.size() == 0) return initialState_;
    else return states_[states_.size() - 1];
  }
};

/**
 * @brief Interface for simulations.
 *
 * A mutation process defines the rules for mutations to occur.
 * The MutationProcess interface provides two methods, one for mutating a character in
 * state i in another character, another for achieving this task n times.
 */
class MutationProcess
{
public:
  MutationProcess() {}
  virtual ~MutationProcess() {}

public:
  /**
   * @brief Mutate a character in state i.
   *
   * @param state The current state of the character.
   */
  virtual size_t mutate(size_t state) const = 0;

  /**
   * @brief Mutate a character in state i n times.
   *
   * @param state The current state of the character.
   * @param n The number of mutations to perform.
   */
  virtual size_t mutate(size_t state, unsigned int n) const = 0;

  /**
   * @brief Get the time before next mutation event.
   *
   * @param state The actual state of the chain;
   * @return A random time before next mutation event.
   */
  virtual double getTimeBeforeNextMutationEvent(size_t state) const = 0;

  /**
   * @brief Simulation a character evolution during a specified time
   * according to the given substitution model and send the final state.
   *
   * @param initialState The state before beginning evolution.
   * @param time         The time during which evolution must occur.
   * @return The resulting state after evolution is completed.
   */
  virtual size_t evolve(size_t initialState, double time) const = 0;

  /**
   * @brief Simulation a character evolution during a specified time
   * according to the given substitution model and send the total path
   * with all intermediate states and times between mutation events.
   *
   * @param initialState The state before beginning evolution.
   * @param time         The time during which evolution must occur.
   * @return The resulting mutation path.
   */
  virtual MutationPath detailedEvolve(size_t initialState, double time) const = 0;

  /**
   * @brief Simulation a character evolution during a specified time
   * according to the given substitution model and send the total
   * path with all intermediate states and times between mutation
   * events, conditional to the final state.q
   *
   * @param initialState The state before beginning evolution.
   * @param finalState The state after  evolution.
   * @param time         The time during which evolution must occur.
   * @return The resulting mutation path.
   */
  virtual MutationPath detailedEvolve(size_t initialState, size_t finalState, double time) const = 0;

  /**
   * @brief Get the substitution model associated to the mutation process.
   *
   * @return The SubstitutionModel associated to this instance.
   */
  virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const = 0;
};

/**
 * @brief Partial implementation of the MutationProcess interface.
 *
 * This class provides an implementation of the MutationProcess interface.
 * It assumes that there are size_ states allowed for the character of interest,
 * and that the distribution of probabilities are in repartition_.
 * As a matter of facts, probabilities must be cumulative, so that repartition_
 * contains values of the repartition function.
 * The mutate function hence draws a random number between 0 and 1 and gives the
 * corresponding character using the bijection of the repartition function.
 *
 * All derived classes must initialize the repartition_ and size_ fields.
 */
class AbstractMutationProcess :
  public virtual MutationProcess
{
protected:
  /**
   * @brief The substitution model to use:
   */
  std::shared_ptr<const SubstitutionModelInterface> model_;

  /**
   * @brief The number of states allowed for the character to mutate.
   */
  size_t size_;

  /**
   * @brief The repartition function for states probabilities.
   *
   * repartition_[i][j] = probability that, being in state i at time t,
   * we'll be in state <= j at time t+1.
   */
  VVdouble repartition_;

public:
  AbstractMutationProcess(std::shared_ptr<const SubstitutionModelInterface> model) :
    model_(model), size_(), repartition_()
  {}

public:
  size_t mutate(size_t state) const;
  size_t mutate(size_t state, unsigned int n) const;
  double getTimeBeforeNextMutationEvent(size_t state) const;
  size_t evolve(size_t initialState, double time) const;
  MutationPath detailedEvolve(size_t initialState, double time) const;
  MutationPath detailedEvolve(size_t initialState, size_t finalState, double time) const;
  std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const { return model_; }
};

/**
 * @brief Generally used mutation process model.
 *
 * This builds a MutationProcess according to a given SubstitutionModel.
 * The underlying mutation process is the following:
 * <ol>
 * <li>Draw a random time @f$ t @f$ from an exponential law with parameter
 * @f$ - \lambda_i @f$,</li>
 * <li> Mutate the initial state. The probability of mutating state @f$i@f$
 * to state @f$j@f$ is:
 * @f[ \frac{Q_{i,j}}{\sum_k Q_{i,k}}. @f]</li>
 * </ol>
 */
class SimpleMutationProcess :
  public AbstractMutationProcess
{
public:
  // Constructor and destructor:

  /**
   * @brief Build a new SimpleMutationProcess object.
   *
   * @param model The substitution model to use.
   */
  SimpleMutationProcess(std::shared_ptr<const SubstitutionModelInterface> model);

  virtual ~SimpleMutationProcess();

  /**
   * @brief Method redefinition for better performance.
   *
   * @param initialState The state before beginning evolution.
   * @param time         The time during which evolution must occur.
   * @return The resulting state after evolution is completed.
   */
  size_t evolve(size_t initialState, double time) const;
};

/**
 * @brief This class is mainly for testing purpose.
 * It allow "self" mutation of the kind i->i;
 */
class SelfMutationProcess : public AbstractMutationProcess
{
public:
  SelfMutationProcess(size_t alphabetSize);

  virtual ~SelfMutationProcess();
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_MUTATIONPROCESS_H
