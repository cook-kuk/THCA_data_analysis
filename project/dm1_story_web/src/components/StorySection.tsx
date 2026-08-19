import React from "react";

interface Props {
  id: string;
  num: string;
  title: string;
  lead?: string;
  children: React.ReactNode;
}

export const StorySection: React.FC<Props> = ({ id, num, title, lead, children }) => (
  <section className="story" id={id}>
    <div className="num">{num}</div>
    <h2>{title}</h2>
    {lead && <div className="lead">{lead}</div>}
    <div className="body">{children}</div>
  </section>
);
