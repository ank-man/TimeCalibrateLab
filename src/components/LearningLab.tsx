/**
 * LearningLab — Structured interactive lesson system.
 * Turns TimeCalibrateLab into a guided educational experience with
 * modules, lessons, challenges, and knowledge checks.
 */

import { useState } from 'react';
import { useStore } from '../store/useStore';
import {
  ChevronRight, ChevronLeft, Check,
  FlaskConical, Lightbulb, Target, ArrowRight, X,
} from 'lucide-react';

export interface LessonStep {
  title: string;
  content: string;
  action?: string;
  hint?: string;
  quiz?: { question: string; options: string[]; correct: number; explanation: string };
}

export interface Lesson {
  id: string;
  title: string;
  description: string;
  icon: string;
  steps: LessonStep[];
  onStart?: () => void;
}

export interface Module {
  id: string;
  title: string;
  description: string;
  color: string;
  lessons: Lesson[];
}

function buildModules(store: ReturnType<typeof useStore.getState>): Module[] {
  return [
    {
      id: 'foundations',
      title: 'Module 1: Foundations',
      description: 'Understand why we need time calibration and the core Bayesian equation.',
      color: 'blue',
      lessons: [
        {
          id: 'why-time',
          title: 'Why Time-Calibrate?',
          description: 'From branch lengths to geological time',
          icon: '🕰️',
          steps: [
            {
              title: 'The Problem',
              content: 'Look at the tree on the left. It shows a primate phylogeny with 8 species. The branch lengths represent **substitutions per site** — how much DNA has changed.\n\nBut substitutions per site is not *time*. A branch of length 0.05 could represent 5 million years (at a rate of 0.01 subst/site/Myr) or 50 million years (at a rate of 0.001).\n\n**We need additional information to convert molecular distances into absolute time.**',
            },
            {
              title: 'The Three Ingredients',
              content: 'To date a tree, we need:\n\n**1. Sequence data** — provides the likelihood (how much evolution happened)\n\n**2. A clock model** — specifies the relationship between substitutions and time (d = rate × time)\n\n**3. Fossil calibrations** — anchor points that connect molecular evolution to geological time\n\nBayes\' theorem combines them:\n\n`Posterior ∝ Likelihood × Prior`',
            },
            {
              title: 'The Bayesian Equation',
              content: 'More precisely:\n\n`P(times, rates | data, fossils) ∝ P(data | times, rates) × P(times | fossils) × P(rates)`\n\n- **P(data | times, rates)** = Likelihood — what the sequences tell us\n- **P(times | fossils)** = Prior on node ages — what fossils tell us\n- **P(rates)** = Prior on evolutionary rates\n- **P(times, rates | data, fossils)** = Posterior — what we want\n\nEverything in this app demonstrates this equation.',
              quiz: {
                question: 'In Bayesian time calibration, what does the POSTERIOR represent?',
                options: [
                  'The probability of the data given our model',
                  'Our belief about node ages BEFORE seeing data',
                  'Our updated belief about node ages AFTER combining data and priors',
                  'The substitution rate',
                ],
                correct: 2,
                explanation: 'The posterior is our updated belief after combining the likelihood (from sequence data) with the prior (from fossils and clock model). It\'s what we actually want to estimate.',
              },
            },
          ],
        },
        {
          id: 'explore-tree',
          title: 'Exploring the Tree',
          description: 'Click nodes, read tooltips, understand the topology',
          icon: '🌳',
          steps: [
            {
              title: 'The Primate Tree',
              content: 'The tree shows 8 primate species: Human, Chimp, Gorilla, Orangutan, Gibbon, Old World Monkey, New World Monkey, and Lemur.\n\nEach **tip** (leaf) represents a living species. Each **internal node** represents a common ancestor — a species that split into two lineages.',
              action: 'Click on different nodes in the tree. Notice which ones have a 🦴 fossil icon.',
            },
            {
              title: 'Node Ages',
              content: 'When you hover over an internal node, you see its **age** in millions of years (Ma).\n\nThe Human-Chimp ancestor (N1) is ~7 Ma. The root of all primates (N7) is ~65 Ma.\n\nThese are the parameters we want to estimate. Right now they show the "true" ages — after running MCMC, they\'ll show posterior estimates with uncertainty.',
              action: 'Click on the node labeled "Human-Chimp" (N1) in the tree. You should see it highlighted.',
            },
            {
              title: 'Calibrated Nodes',
              content: 'Three nodes have fossil calibrations (🦴 icon):\n\n- **Human-Chimp (N1)** — Lognormal calibration based on Sahelanthropus (~6-7 Ma)\n- **Catarrhini (N5)** — Uniform calibration (25-35 Ma)\n- **Primates root (N7)** — Exponential calibration (offset 55 Ma)\n\nThese calibrations are our **priors** — what we believe about node ages before looking at sequence data.',
              quiz: {
                question: 'Why do calibrated nodes have different distribution types?',
                options: [
                  'It doesn\'t matter — all distributions give the same answer',
                  'Different fossil evidence supports different shapes of uncertainty',
                  'MCMCTree requires different types for different nodes',
                  'Lognormal is always better',
                ],
                correct: 1,
                explanation: 'The distribution shape should reflect the nature of your fossil evidence. A fossil gives a firm minimum age (the clade is at least this old) but the maximum is uncertain — lognormal captures this asymmetry. When you only know a range, uniform is appropriate.',
              },
            },
          ],
        },
      ],
    },
    {
      id: 'priors',
      title: 'Module 2: Priors & Calibrations',
      description: 'Learn how fossil calibrations shape your prior beliefs about node ages.',
      color: 'amber',
      lessons: [
        {
          id: 'prior-types',
          title: 'Calibration Distributions',
          description: 'Uniform, Lognormal, and Exponential — when to use each',
          icon: '📊',
          steps: [
            {
              title: 'Uniform Distribution',
              content: '**Uniform(min, max)** assigns equal probability to all ages within the range.\n\nUse when: You know the minimum (fossil) and maximum (geological event) but have no reason to prefer any age within the range.\n\nExample: "This clade is between 25 and 35 Ma old" → Uniform(25, 35)',
              action: 'In the Calibration tab on the right, select node "Catarrhini (N5)" and note it has a Uniform distribution. Try adjusting the min and max sliders.',
              hint: 'Watch the density plot update as you move the sliders.',
            },
            {
              title: 'Lognormal Distribution',
              content: '**Lognormal(offset, μ, σ)** is right-skewed. Most probability is near the offset (fossil age), with a long tail toward older ages.\n\nUse when: You have a good fossil minimum but the true age could be substantially older.\n\nThe offset is the fossil age (hard minimum). μ and σ control the shape above the offset.',
              action: 'Select "Human-Chimp (N1)" which has a Lognormal calibration. Increase σ to see the distribution become wider and more uncertain.',
            },
            {
              title: 'Exponential Distribution',
              content: '**Exponential(offset, λ)** puts most probability mass right at the offset and decays quickly.\n\nUse when: You believe the true age is very close to the fossil age — the fossil record for this group is good.\n\nSmaller λ → wider distribution → more uncertainty above the fossil age.',
              action: 'Select "Primates (N7)" which has an Exponential calibration. Try different λ values and see how the density changes.',
              quiz: {
                question: 'If you have a fossil at 65 Ma but think the clade could easily be 80-100 Ma old, which distribution is most appropriate?',
                options: [
                  'Exponential(65, 1.0) — tight around the fossil',
                  'Uniform(65, 65) — exactly 65 Ma',
                  'Lognormal(65, 1.5, 0.8) — peaks above fossil, long tail',
                  'Uniform(0, 200) — we have no idea',
                ],
                correct: 2,
                explanation: 'A lognormal with offset at 65 Ma would place the mode above the fossil age with a right tail extending toward older ages. This captures the idea that 65 Ma is a minimum but the true age is likely somewhat older.',
              },
            },
            {
              title: 'The Width Matters',
              content: '**Wide calibrations** (e.g., Uniform 50-150 Ma) express high uncertainty. The prior is spread over a large range.\n\n**Narrow calibrations** (e.g., Uniform 63-67 Ma) express confidence. The prior is concentrated.\n\n**Key insight:** When your sequence data is weak (short alignment), the posterior will look like the prior. When data is strong (long alignment), the posterior will be driven by the likelihood regardless of the prior width.',
              action: 'Try this: Go to the Scenarios menu (top bar) and load "Wide Calibration". Then click "Run MCMC". Notice how the posterior is dominated by the prior.',
            },
          ],
          onStart: () => {
            store.selectNode('N5');
          },
        },
        {
          id: 'prior-conflicts',
          title: 'Conflicting Calibrations',
          description: 'What happens when calibrations contradict each other',
          icon: '⚠️',
          steps: [
            {
              title: 'How Conflicts Arise',
              content: 'A **calibration conflict** occurs when the prior for a child node overlaps with or exceeds the prior for its parent node.\n\nExample: If the parent is calibrated as Uniform(18, 22) and the child as Uniform(20, 25), then the child could be *older* than the parent — which is biologically impossible.',
            },
            {
              title: 'See It In Action',
              content: 'Let\'s load a scenario with conflicting calibrations and observe the warning.',
              action: 'Go to Scenarios → "Conflicting Calibrations". Notice the red warning banner that appears in the Calibration panel.',
            },
            {
              title: 'What MCMCTree Does',
              content: 'When calibrations conflict, the MCMC sampler has trouble:\n\n- Many proposals are rejected (low acceptance rate)\n- The posterior is squeezed into a narrow sliver where both constraints can be satisfied\n- Results may depend heavily on which constraint "wins"\n\n**In real analyses:** Always check for conflicts before running MCMCTree. Run a prior-only analysis (usedata=0) and verify the effective prior makes biological sense.',
              quiz: {
                question: 'You get an acceptance rate of 3% in your MCMCTree run. What is the most likely cause?',
                options: [
                  'Your alignment is too long',
                  'Conflicting calibrations or impossible age combinations',
                  'The strict clock is correct',
                  'You need more MCMC iterations',
                ],
                correct: 1,
                explanation: 'Very low acceptance rates usually indicate that most proposed node ages violate constraints — typically from conflicting calibrations or overly tight bounds. The sampler keeps proposing impossible states and rejecting them.',
              },
            },
          ],
          onStart: () => {
            store.loadScenario('conflicting_calibrations');
          },
        },
      ],
    },
    {
      id: 'clock-models',
      title: 'Module 3: The Molecular Clock',
      description: 'Strict vs relaxed clocks and how rate variation affects inference.',
      color: 'purple',
      lessons: [
        {
          id: 'strict-vs-relaxed',
          title: 'Strict vs Relaxed Clock',
          description: 'See how clock choice changes everything',
          icon: '⏱️',
          steps: [
            {
              title: 'The Strict Clock',
              content: '**Strict clock:** Every branch evolves at exactly the same rate.\n\n`distance = rate × time` with one global rate for all branches.\n\nThis was the original molecular clock hypothesis (Zuckerkandl & Pauling, 1962). It\'s elegant but almost never true — mice evolve faster than elephants, parasites faster than hosts.',
              action: 'Click the "Clock" tab on the right panel. Select "Strict Clock". Then Run MCMC.',
            },
            {
              title: 'The Relaxed Clock',
              content: '**Relaxed clock:** Each branch gets its own rate drawn from a distribution:\n\n`r_i ~ LogNormal(μ, σ²)`\n\nThe parameter **σ** controls how much rates vary:\n- σ ≈ 0 → essentially strict clock\n- σ = 0.1 → mild variation\n- σ = 1.0 → extreme variation\n\nHigher σ means more rate heterogeneity, which means more uncertainty in divergence time estimates.',
              action: 'Switch to "Relaxed Clock". Notice the σ slider appears. Run MCMC and compare the credible intervals to the strict clock results.',
            },
            {
              title: 'The σ Effect',
              content: 'This is the single most important parameter for the width of your credible intervals.\n\n**Low σ** → branches have similar rates → times are well-constrained\n**High σ** → branches have wildly different rates → times are uncertain\n\nIn MCMCTree, σ² is controlled by the `sigma2_gamma` prior.',
              action: 'Load the "High Relaxed Variance" scenario from the menu. Run MCMC. See how the posteriors become much wider compared to the strict clock.',
              quiz: {
                question: 'You switch from strict to relaxed clock and your 95% HPD intervals double in width. Why?',
                options: [
                  'The relaxed clock is wrong for your data',
                  'Rate uncertainty propagates to time uncertainty — if rates vary, the same branch length maps to a range of times',
                  'You need more MCMC iterations',
                  'The alignment is too short',
                ],
                correct: 1,
                explanation: 'Under a relaxed clock, the same observed branch length (substitutions) could have been produced by a fast rate over a short time OR a slow rate over a long time. This rate-time ambiguity widens the credible intervals.',
              },
            },
            {
              title: 'Branch Rate Heatmap',
              content: 'When using the relaxed clock, after running MCMC, look at the **Branch Rate Heatmap** in the Clock panel.\n\nIt shows the estimated rate for each branch. Blue = slow, red = fast.\n\nIn real data, you might see:\n- Terminal branches to rodents are fast (short generation time)\n- Deep branches are intermediate\n- Some lineages are conspicuously slow (e.g., coelacanth)',
              action: 'With relaxed clock results visible, examine the Branch Rate Heatmap. Which lineages have the fastest/slowest rates?',
            },
          ],
          onStart: () => {
            store.setClockType('strict');
          },
        },
      ],
    },
    {
      id: 'data-power',
      title: 'Module 4: The Power of Data',
      description: 'How alignment length and substitution model affect inference.',
      color: 'green',
      lessons: [
        {
          id: 'alignment-length',
          title: 'Alignment Length Effect',
          description: 'See how more data tightens the posterior',
          icon: '🧬',
          steps: [
            {
              title: 'Weak Data: 200 Sites',
              content: 'With a short alignment (200 sites), the likelihood is **weak**. It provides a broad, flat signal about divergence times.\n\nResult: The posterior is dominated by the prior. Your fossil calibrations determine the answer.',
              action: 'Go to the Data tab. Set alignment length to 200. Run MCMC. Look at the density plot for any calibrated node — the posterior (blue) closely follows the prior (amber dashed).',
            },
            {
              title: 'Strong Data: 10,000 Sites',
              content: 'With a long alignment (10,000 sites), the likelihood is **strong**. It provides a sharp, peaked signal.\n\nResult: The posterior is driven by the data. Even wide priors get overridden.',
              action: 'Now increase alignment length to 10,000. Run MCMC again. The posterior should be much narrower than the prior.',
              hint: 'This is Bayes\' theorem in action: strong data overwhelms the prior.',
            },
            {
              title: 'The Crossover',
              content: 'Somewhere between 200 and 10,000 sites, there\'s a crossover point where data and priors contribute equally.\n\n**This is why genomic-scale datasets changed everything.** With thousands of genes and millions of sites, sequence data dominates. But for poorly supported nodes (low signal), calibrations still matter.',
              quiz: {
                question: 'You have Uniform(50, 150) as your root calibration. With 100 sites, the posterior mean is 98 Ma. With 10,000 sites, it\'s 72 Ma. What happened?',
                options: [
                  'The 10,000-site result is wrong',
                  'With more data, the likelihood peak at ~72 Ma pulled the posterior away from the prior mean (~100 Ma)',
                  'The calibration changed',
                  'MCMC didn\'t converge in the first run',
                ],
                correct: 1,
                explanation: 'With 100 sites, the likelihood was flat and the posterior tracked the prior (mean ≈ center of Uniform = 100). With 10,000 sites, the likelihood was sharp and peaked at ~72 Ma, pulling the posterior toward the data-driven estimate.',
              },
            },
          ],
          onStart: () => {
            store.loadScenario('wide_calibration');
          },
        },
        {
          id: 'likelihood-models',
          title: 'JC69 vs Distance Likelihood',
          description: 'Compare substitution models',
          icon: '📐',
          steps: [
            {
              title: 'JC69 Model',
              content: 'The **JC69 model** (Jukes-Cantor 1969) is the simplest nucleotide substitution model:\n\n- All bases equally frequent (25% each)\n- All substitution types equally likely\n- P(same base after distance d) = 1/4 + 3/4 × e^(-4d/3)\n\nIt\'s simple but captures the essential feature: the probability of observing a difference increases with evolutionary distance but saturates at 75%.',
            },
            {
              title: 'Distance Approximation',
              content: 'The **distance approximation** uses a Gaussian likelihood:\n\n`observed_distance ~ Normal(expected_distance, variance/L)`\n\nThis is faster to compute and works well when distances are not saturated. It\'s essentially a Taylor expansion of the JC69 likelihood.',
              action: 'In the Data tab, toggle between JC69 and Distance mode. Run MCMC with each. The results should be very similar for moderate divergences.',
            },
            {
              title: 'When Models Matter',
              content: 'For divergence time estimation, the substitution model choice (JC69 vs HKY vs GTR) has **less effect** than you might think — the calibrations and clock model typically dominate.\n\nHowever, for very divergent sequences (>0.5 subst/site), model choice matters more because saturation correction differs between models.\n\n**In MCMCTree:** Use `model = 4` (HKY) or `model = 7` (GTR) for real analyses.',
            },
          ],
        },
      ],
    },
    {
      id: 'mcmc',
      title: 'Module 5: MCMC Sampling',
      description: 'Understand how Markov Chain Monte Carlo explores the posterior.',
      color: 'rose',
      lessons: [
        {
          id: 'mcmc-basics',
          title: 'How MCMC Works',
          description: 'Proposals, acceptance, convergence',
          icon: '🎲',
          steps: [
            {
              title: 'The Idea',
              content: 'We want to know the posterior distribution but can\'t compute it directly. MCMC generates **samples** from the posterior without knowing its exact shape.\n\n**Metropolis-Hastings algorithm:**\n1. Start with initial node ages\n2. Propose a small random change to one age\n3. Compute acceptance ratio: α = P(new)/P(old)\n4. Accept with probability α, else keep old values\n5. Record the current state\n6. Repeat thousands of times',
            },
            {
              title: 'Run and Watch',
              content: 'When you click "Run MCMC", the app runs 3,000 iterations of Metropolis-Hastings in a Web Worker (off the main thread so the UI stays responsive).\n\nAfter burn-in (first 25% discarded), the remaining samples approximate the posterior.',
              action: 'Click "Run MCMC" and watch the progress bar. When done, check the acceptance rate in the summary bar below the toolbar.',
            },
            {
              title: 'Acceptance Rate',
              content: 'The **acceptance rate** tells you how often proposals were accepted:\n\n- **20-40%** → healthy (optimal ~23%)\n- **< 10%** → proposals are too large, most get rejected\n- **> 70%** → proposals are too small, chain explores slowly\n\nIn MCMCTree, the `finetune` parameter controls proposal step sizes, and `finetune = 1:` enables auto-tuning.',
            },
            {
              title: 'Trace Plots',
              content: 'Switch to **Advanced mode** (top bar) to see the **log-posterior trace plot** at the bottom of the right panel.\n\nA good trace looks like a "hairy caterpillar" — fluctuating around a stable value. If it\'s trending upward, the chain hasn\'t converged yet.',
              action: 'Switch to Advanced mode. Run MCMC. Look at the trace plot. Does it look stationary?',
              quiz: {
                question: 'Your trace plot shows the log-posterior climbing for the first 50,000 iterations then stabilizing. What should you do?',
                options: [
                  'The run is fine, use all samples',
                  'Set burnin > 50,000 to discard the non-stationary part',
                  'Reduce the number of samples',
                  'Switch to strict clock',
                ],
                correct: 1,
                explanation: 'The climbing phase is the chain finding the high-probability region. Those samples don\'t represent the posterior. Set burnin to at least 50,000 (preferably more) to discard them. Only the stable, stationary samples should be used for inference.',
              },
            },
          ],
        },
        {
          id: 'hpd-intervals',
          title: 'HPD Intervals',
          description: 'Understanding highest posterior density intervals',
          icon: '📏',
          steps: [
            {
              title: 'What is an HPD?',
              content: 'The **95% Highest Posterior Density (HPD) interval** is the shortest interval containing 95% of the posterior probability.\n\nUnlike a symmetric 95% CI, the HPD is the *shortest* possible interval — it always includes the highest-density (most probable) region.\n\nFor skewed posteriors, the HPD is asymmetric.',
            },
            {
              title: 'Reading HPD Intervals',
              content: 'After running MCMC, hover over any internal node. You\'ll see:\n\n- **Mean age** — the average of all posterior samples\n- **95% HPD** — the range within which the true age lies with 95% probability\n\nWide HPDs mean high uncertainty. Narrow HPDs mean the data and priors strongly constrain the age.',
              action: 'Run MCMC, then click different internal nodes. Compare HPD widths between calibrated and uncalibrated nodes.',
            },
            {
              title: 'What Widens HPDs?',
              content: '**Factors that widen HPDs:**\n- Short alignment (weak likelihood)\n- Wide calibration priors\n- High rate variation (relaxed clock with large σ)\n- Few calibration points\n- Node far from calibrated nodes\n\n**Factors that narrow HPDs:**\n- Long alignment (strong likelihood)\n- Tight calibration priors\n- Strict or near-strict clock\n- Multiple nearby calibrations',
              quiz: {
                question: 'Node X has HPD [25, 95] Ma. What does this NOT mean?',
                options: [
                  'There is 95% posterior probability that the true age is between 25 and 95 Ma',
                  'The true age is definitely between 25 and 95 Ma',
                  'The data and priors are consistent with ages in this range',
                  'This is the shortest interval containing 95% of the posterior mass',
                ],
                correct: 1,
                explanation: 'The HPD is a credible interval — it says there\'s 95% posterior probability within this range, but there\'s still 5% probability outside it. It does NOT guarantee the true age is inside. Also, if your model or calibrations are wrong, the interval may not contain the truth at all.',
              },
            },
          ],
        },
      ],
    },
    {
      id: 'putting-together',
      title: 'Module 6: Putting It All Together',
      description: 'Combine everything to understand how real analyses work.',
      color: 'teal',
      lessons: [
        {
          id: 'full-workflow',
          title: 'A Complete Analysis',
          description: 'Walk through the full Bayesian dating workflow',
          icon: '🔬',
          steps: [
            {
              title: 'Step 1: Set Calibrations',
              content: 'In real analyses, you start by reviewing the fossil record for your clade and selecting appropriate calibration points.\n\n**Best practices:**\n- Use 3-5 well-justified calibrations spread across the tree\n- Prefer calibrations with clear phylogenetic placement\n- Use the distribution that matches your fossil evidence uncertainty\n- Check for conflicts between calibrations',
              action: 'Reset the app. Examine the default calibrations on N1, N5, and N7. These represent typical primate calibrations.',
            },
            {
              title: 'Step 2: Choose Clock Model',
              content: 'For most real datasets, use the **relaxed clock** (independent rates, MCMCTree clock=2).\n\nAlso run with autocorrelated rates (clock=3) as a sensitivity check.\n\n**Rule of thumb:** If strict and relaxed give very different answers, the strict clock is probably inappropriate.',
            },
            {
              title: 'Step 3: Run Prior-Only',
              content: 'Before using data, run with `usedata = 0` to see the **effective prior** — what the software actually uses after combining your calibrations with the birth-death tree prior.\n\nIf the effective prior differs substantially from your intended calibrations, adjust BDparas or calibrations.',
            },
            {
              title: 'Step 4: Run with Data',
              content: 'Run the full analysis with `usedata = 2` (approximate likelihood). Run at least 2 independent chains with different random seeds.\n\nCheck convergence:\n- ESS > 200 for all parameters\n- Trace plots stationary\n- Independent chains agree',
            },
            {
              title: 'Step 5: Interpret',
              content: 'Report **posterior means with 95% HPD intervals**.\n\nWide intervals are honest — they reflect genuine uncertainty. Don\'t be tempted to use a tighter prior just to get narrower intervals.\n\nCompare results under different models and calibration strategies. Robust conclusions hold across analyses.',
              quiz: {
                question: 'You run MCMCTree with two independent chains. Chain 1 gives root age 72 ± 15 Ma, Chain 2 gives 91 ± 20 Ma. What do you conclude?',
                options: [
                  'Average them: root age is ~82 Ma',
                  'The chains have NOT converged — run much longer before trusting results',
                  'Use whichever chain has higher ESS',
                  'The true age is between 72 and 91 Ma',
                ],
                correct: 1,
                explanation: 'If independent chains give substantially different posteriors, the MCMC has NOT converged. No result from either chain is trustworthy. You need to run much longer (10x or more) or diagnose why convergence failed (conflicting calibrations, poor proposal tuning, etc).',
              },
            },
          ],
          onStart: () => {
            store.initialize();
          },
        },
      ],
    },
  ];
}

export default function LearningLab({ onClose }: { onClose: () => void }) {
  const store = useStore();
  const theme = store.theme;
  const modules = buildModules(useStore.getState());

  const [activeModuleIdx, setActiveModuleIdx] = useState<number | null>(null);
  const [activeLessonIdx, setActiveLessonIdx] = useState(0);
  const [activeStepIdx, setActiveStepIdx] = useState(0);
  const [completedSteps, setCompletedSteps] = useState<Set<string>>(new Set());
  const [quizAnswer, setQuizAnswer] = useState<number | null>(null);
  const [showQuizResult, setShowQuizResult] = useState(false);

  const isDark = theme === 'dark';
  const bg = isDark ? 'bg-slate-900' : 'bg-gray-50';
  const cardBg = isDark ? 'bg-slate-800' : 'bg-white';
  const border = isDark ? 'border-slate-700' : 'border-gray-200';
  const textMuted = isDark ? 'text-slate-400' : 'text-gray-500';

  const activeModule = activeModuleIdx !== null ? modules[activeModuleIdx] : null;
  const activeLesson = activeModule ? activeModule.lessons[activeLessonIdx] : null;
  const activeStep = activeLesson ? activeLesson.steps[activeStepIdx] : null;

  function stepKey(mi: number, li: number, si: number) {
    return `${mi}-${li}-${si}`;
  }

  function markComplete() {
    if (activeModuleIdx === null) return;
    const key = stepKey(activeModuleIdx, activeLessonIdx, activeStepIdx);
    setCompletedSteps((prev) => new Set(prev).add(key));
  }

  function nextStep() {
    if (!activeLesson || !activeModule || activeModuleIdx === null) return;
    markComplete();
    setQuizAnswer(null);
    setShowQuizResult(false);

    if (activeStepIdx < activeLesson.steps.length - 1) {
      setActiveStepIdx(activeStepIdx + 1);
    } else if (activeLessonIdx < activeModule.lessons.length - 1) {
      setActiveLessonIdx(activeLessonIdx + 1);
      setActiveStepIdx(0);
      const nextLesson = activeModule.lessons[activeLessonIdx + 1];
      if (nextLesson.onStart) nextLesson.onStart();
    } else if (activeModuleIdx < modules.length - 1) {
      setActiveModuleIdx(activeModuleIdx + 1);
      setActiveLessonIdx(0);
      setActiveStepIdx(0);
    }
  }

  function prevStep() {
    setQuizAnswer(null);
    setShowQuizResult(false);
    if (activeStepIdx > 0) {
      setActiveStepIdx(activeStepIdx - 1);
    } else if (activeLessonIdx > 0) {
      setActiveLessonIdx(activeLessonIdx - 1);
      const prevLesson = activeModule!.lessons[activeLessonIdx - 1];
      setActiveStepIdx(prevLesson.steps.length - 1);
    }
  }

  function startModule(idx: number) {
    setActiveModuleIdx(idx);
    setActiveLessonIdx(0);
    setActiveStepIdx(0);
    setQuizAnswer(null);
    setShowQuizResult(false);
    const firstLesson = modules[idx].lessons[0];
    if (firstLesson.onStart) firstLesson.onStart();
  }

  function goToOverview() {
    setActiveModuleIdx(null);
    setActiveLessonIdx(0);
    setActiveStepIdx(0);
    setQuizAnswer(null);
    setShowQuizResult(false);
  }

  const totalSteps = modules.reduce((sum, m) => sum + m.lessons.reduce((s, l) => s + l.steps.length, 0), 0);
  const completedCount = completedSteps.size;
  const progressPct = totalSteps > 0 ? (completedCount / totalSteps) * 100 : 0;

  // Module overview
  if (activeModuleIdx === null) {
    return (
      <div className={`fixed inset-0 z-50 flex flex-col ${bg} overflow-y-auto`}>
        <div className={`sticky top-0 z-10 flex items-center justify-between px-6 py-3 border-b ${border} ${cardBg} backdrop-blur-sm`}>
          <div className="flex items-center gap-3">
            <FlaskConical size={22} className="text-blue-500" />
            <div>
              <h1 className="text-lg font-bold">TimeCalibrateLab Learning Lab</h1>
              <p className={`text-xs ${textMuted}`}>Interactive lessons on Bayesian phylogenetic time calibration</p>
            </div>
          </div>
          <div className="flex items-center gap-4">
            <div className="text-xs text-right">
              <div className={textMuted}>Progress</div>
              <div className="font-semibold">{completedCount}/{totalSteps} steps</div>
            </div>
            <div className="w-32 h-2 rounded-full bg-slate-700 overflow-hidden">
              <div className="h-full bg-blue-500 rounded-full transition-all" style={{ width: `${progressPct}%` }} />
            </div>
            <button onClick={onClose} className={`p-2 rounded-lg hover:bg-slate-700 ${textMuted}`} title="Exit Learning Lab">
              <X size={18} />
            </button>
          </div>
        </div>

        <div className="max-w-4xl mx-auto w-full px-6 py-8">
          <div className="text-center mb-10">
            <h2 className="text-2xl font-bold mb-2">Learn Bayesian Time Calibration</h2>
            <p className={`${textMuted} max-w-lg mx-auto`}>
              Six interactive modules that take you from zero to understanding every aspect of
              molecular clock dating. Each lesson includes hands-on exercises with the simulator.
            </p>
          </div>

          <div className="grid gap-4">
            {modules.map((mod, mi) => {
              const modSteps = mod.lessons.reduce((s, l) => s + l.steps.length, 0);
              const modCompleted = mod.lessons.reduce((s, l, li) =>
                s + l.steps.filter((_, si) => completedSteps.has(stepKey(mi, li, si))).length, 0);
              const colorMap: Record<string, string> = {
                blue: 'from-blue-600 to-blue-700',
                amber: 'from-amber-600 to-amber-700',
                purple: 'from-purple-600 to-purple-700',
                green: 'from-green-600 to-green-700',
                rose: 'from-rose-600 to-rose-700',
                teal: 'from-teal-600 to-teal-700',
              };
              const grad = colorMap[mod.color] || colorMap.blue;
              const isComplete = modCompleted === modSteps && modSteps > 0;

              return (
                <button
                  key={mod.id}
                  onClick={() => startModule(mi)}
                  className={`w-full text-left p-5 rounded-xl border ${border} ${cardBg} hover:shadow-lg transition-all group`}
                  title={`Start ${mod.title}`}
                >
                  <div className="flex items-start gap-4">
                    <div className={`w-12 h-12 rounded-xl bg-gradient-to-br ${grad} flex items-center justify-center text-white font-bold text-sm shrink-0`}>
                      {isComplete ? <Check size={20} /> : `M${mi + 1}`}
                    </div>
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2">
                        <h3 className="font-bold">{mod.title}</h3>
                        <span className={`text-xs ${textMuted}`}>{mod.lessons.length} lessons</span>
                      </div>
                      <p className={`text-sm ${textMuted} mt-0.5`}>{mod.description}</p>
                      <div className="flex items-center gap-2 mt-2">
                        {mod.lessons.map((l) => (
                          <span key={l.id} className={`text-xs px-2 py-0.5 rounded-full ${isDark ? 'bg-slate-700' : 'bg-gray-100'}`}>
                            {l.icon} {l.title}
                          </span>
                        ))}
                      </div>
                    </div>
                    <div className="flex items-center gap-2 shrink-0">
                      <span className={`text-xs ${textMuted}`}>{modCompleted}/{modSteps}</span>
                      <ArrowRight size={16} className={`${textMuted} group-hover:translate-x-1 transition-transform`} />
                    </div>
                  </div>
                </button>
              );
            })}
          </div>
        </div>
      </div>
    );
  }

  // Lesson view
  return (
    <div className={`fixed inset-y-0 right-0 w-[420px] z-50 flex flex-col ${cardBg} border-l ${border} shadow-2xl`}>
      {/* Header */}
      <div className={`flex items-center justify-between px-4 py-2.5 border-b ${border}`}>
        <div className="flex items-center gap-2 min-w-0">
          <button onClick={goToOverview} className={`p-1 rounded hover:bg-slate-700 ${textMuted} shrink-0`} title="Back to modules">
            <ChevronLeft size={16} />
          </button>
          <div className="min-w-0">
            <div className={`text-[10px] uppercase tracking-wider ${textMuted}`}>{activeModule?.title}</div>
            <div className="text-sm font-semibold truncate">{activeLesson?.icon} {activeLesson?.title}</div>
          </div>
        </div>
        <button onClick={onClose} className={`p-1.5 rounded hover:bg-slate-700 ${textMuted}`} title="Close">
          <X size={16} />
        </button>
      </div>

      {/* Progress dots */}
      <div className={`flex items-center gap-1.5 px-4 py-2 border-b ${border}`}>
        {activeLesson?.steps.map((_, si) => {
          const key = stepKey(activeModuleIdx!, activeLessonIdx, si);
          const done = completedSteps.has(key);
          const current = si === activeStepIdx;
          return (
            <button
              key={si}
              onClick={() => { setActiveStepIdx(si); setQuizAnswer(null); setShowQuizResult(false); }}
              className={`w-6 h-1.5 rounded-full transition-all ${
                current ? 'bg-blue-500' : done ? 'bg-green-500' : isDark ? 'bg-slate-700' : 'bg-gray-300'
              }`}
              title={`Step ${si + 1}`}
            />
          );
        })}
        <span className={`text-[10px] ml-auto ${textMuted}`}>
          {activeStepIdx + 1}/{activeLesson?.steps.length}
        </span>
      </div>

      {/* Content */}
      <div className="flex-1 overflow-y-auto px-4 py-4 space-y-4">
        {activeStep && (
          <div className="animate-fade-in">
            <h3 className="text-base font-bold mb-3">{activeStep.title}</h3>

            <div className={`text-sm leading-relaxed space-y-2 ${isDark ? 'text-slate-300' : 'text-gray-600'}`}>
              {activeStep.content.split('\n\n').map((paragraph, i) => (
                <p key={i} dangerouslySetInnerHTML={{
                  __html: paragraph
                    .replace(/\*\*(.+?)\*\*/g, '<strong class="font-semibold text-white">$1</strong>')
                    .replace(/`(.+?)`/g, `<code class="px-1 py-0.5 rounded text-xs font-mono ${isDark ? 'bg-slate-700 text-blue-300' : 'bg-gray-100 text-blue-600'}">$1</code>`)
                    .replace(/\n/g, '<br/>')
                }} />
              ))}
            </div>

            {/* Action prompt */}
            {activeStep.action && (
              <div className={`mt-4 p-3 rounded-lg border ${isDark ? 'bg-blue-900/20 border-blue-800/40' : 'bg-blue-50 border-blue-200'}`}>
                <div className="flex items-start gap-2">
                  <Target size={14} className="text-blue-400 mt-0.5 shrink-0" />
                  <div>
                    <div className="text-xs font-semibold text-blue-400 mb-0.5">Try it</div>
                    <div className={`text-xs ${isDark ? 'text-blue-200' : 'text-blue-700'}`}>{activeStep.action}</div>
                  </div>
                </div>
              </div>
            )}

            {/* Hint */}
            {activeStep.hint && (
              <div className={`mt-2 p-2 rounded-lg ${isDark ? 'bg-amber-900/20 border border-amber-800/40' : 'bg-amber-50 border border-amber-200'}`}>
                <div className="flex items-start gap-2">
                  <Lightbulb size={12} className="text-amber-400 mt-0.5 shrink-0" />
                  <div className={`text-xs ${isDark ? 'text-amber-200' : 'text-amber-700'}`}>{activeStep.hint}</div>
                </div>
              </div>
            )}

            {/* Quiz */}
            {activeStep.quiz && (
              <div className={`mt-4 p-4 rounded-lg border ${isDark ? 'bg-slate-800/80 border-slate-600' : 'bg-gray-50 border-gray-200'}`}>
                <div className="flex items-center gap-2 mb-3">
                  <FlaskConical size={14} className="text-purple-400" />
                  <span className="text-xs font-bold uppercase tracking-wider text-purple-400">Knowledge Check</span>
                </div>
                <p className="text-sm font-medium mb-3">{activeStep.quiz.question}</p>
                <div className="space-y-1.5">
                  {activeStep.quiz.options.map((opt, oi) => {
                    const isSelected = quizAnswer === oi;
                    const isCorrect = oi === activeStep.quiz!.correct;
                    let btnClass = `w-full text-left px-3 py-2 rounded-lg text-xs border transition-all `;
                    if (showQuizResult) {
                      if (isCorrect) btnClass += 'bg-green-900/30 border-green-600 text-green-300';
                      else if (isSelected && !isCorrect) btnClass += 'bg-red-900/30 border-red-600 text-red-300';
                      else btnClass += `${isDark ? 'border-slate-600 text-slate-400' : 'border-gray-300 text-gray-400'}`;
                    } else {
                      btnClass += isSelected
                        ? 'bg-blue-900/30 border-blue-500 text-blue-300'
                        : `${isDark ? 'border-slate-600 hover:border-slate-500' : 'border-gray-300 hover:border-gray-400'}`;
                    }
                    return (
                      <button
                        key={oi}
                        onClick={() => { if (!showQuizResult) setQuizAnswer(oi); }}
                        className={btnClass}
                        title={`Option ${oi + 1}`}
                      >
                        <span className="font-mono mr-2 opacity-50">{String.fromCharCode(65 + oi)}.</span>
                        {opt}
                        {showQuizResult && isCorrect && <Check size={12} className="inline ml-1 text-green-400" />}
                      </button>
                    );
                  })}
                </div>
                {quizAnswer !== null && !showQuizResult && (
                  <button
                    onClick={() => setShowQuizResult(true)}
                    className="mt-3 px-4 py-1.5 bg-purple-600 text-white rounded-lg text-xs font-medium hover:bg-purple-700 transition-colors"
                    title="Check answer"
                  >
                    Check Answer
                  </button>
                )}
                {showQuizResult && (
                  <div className={`mt-3 p-2 rounded-lg text-xs ${
                    quizAnswer === activeStep.quiz!.correct
                      ? isDark ? 'bg-green-900/20 text-green-200' : 'bg-green-50 text-green-700'
                      : isDark ? 'bg-red-900/20 text-red-200' : 'bg-red-50 text-red-700'
                  }`}>
                    <strong>{quizAnswer === activeStep.quiz!.correct ? 'Correct!' : 'Not quite.'}</strong> {activeStep.quiz!.explanation}
                  </div>
                )}
              </div>
            )}
          </div>
        )}
      </div>

      {/* Navigation */}
      <div className={`flex items-center justify-between px-4 py-3 border-t ${border}`}>
        <button
          onClick={prevStep}
          disabled={activeStepIdx === 0 && activeLessonIdx === 0}
          className={`flex items-center gap-1 px-3 py-1.5 rounded-lg text-xs font-medium transition-colors ${
            activeStepIdx === 0 && activeLessonIdx === 0
              ? 'opacity-30 cursor-not-allowed'
              : `${isDark ? 'hover:bg-slate-700' : 'hover:bg-gray-200'}`
          } ${textMuted}`}
          title="Previous step"
        >
          <ChevronLeft size={14} /> Previous
        </button>
        <button
          onClick={nextStep}
          className="flex items-center gap-1 px-4 py-1.5 bg-blue-600 text-white rounded-lg text-xs font-semibold hover:bg-blue-700 transition-colors"
          title="Next step"
        >
          {activeStepIdx === (activeLesson?.steps.length ?? 1) - 1 &&
           activeLessonIdx === (activeModule?.lessons.length ?? 1) - 1
            ? 'Complete Module'
            : 'Next'
          }
          <ChevronRight size={14} />
        </button>
      </div>
    </div>
  );
}
