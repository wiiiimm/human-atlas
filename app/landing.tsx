export default function Landing() {
  return (
    <main className="atlas-landing">
      <header className="landing-brand">Human Atlas <span>3D</span></header>
      <section className="landing-intro" aria-labelledby="landing-title">
        <p className="landing-eyebrow">EXPLORE HUMAN ANATOMY</p>
        <h1 id="landing-title">Understand the body.<br/>One structure at a time.</h1>
        <p>Explore muscles, bones, and organs in 3D. Choose a model to begin.</p>
      </section>
      <nav className="model-cards" aria-label="Anatomy models">
        <a className="model-card female-card" href="/female">
          <span className="model-card-label">FEMALE</span>
          <h2>Female anatomy</h2>
          <p>A study model with female anatomy and estimated body proportions.</p>
          <span className="model-card-action">Explore female model <span aria-hidden="true">↗</span></span>
        </a>
        <a className="model-card male-card" href="/male">
          <span className="model-card-label">MALE</span>
          <h2>Male anatomy</h2>
          <p>Explore the BodyParts3D model, from whole systems to individual structures.</p>
          <span className="model-card-action">Explore male model <span aria-hidden="true">↗</span></span>
        </a>
      </nav>
      <footer className="landing-footer">Rotate the body. Reveal its layers. Find a structure.</footer>
    </main>
  );
}
